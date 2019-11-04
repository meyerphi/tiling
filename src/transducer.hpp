#pragma once

#include <vector>
#include <stack>
#include <limits>
#include <algorithm>

#include "tiles.hpp"

typedef int Index;

struct Edge {
    Index source;
    Index target;
    Color north;
    Color south;
    std::vector<TileIndex> tiles;

    Edge(Index source, Index target, Color north, Color south) :
        source(source), target(target), north(north), south(south)
    {}

    /*
     * Note: equality function does not compare the set of tiles,
     * as it is sufficient to keep one edge with otherwise same values.
     */
    bool operator==(const Edge& rhs) const {
        return source == rhs.source && target == rhs.target &&
               north  == rhs.north  && south  == rhs.south ;
    }
    bool operator!=(const Edge& rhs) const {
        return !operator==(rhs);
    }
};

// optimized stack frame for Tarjan and DFS
struct StackFrame {
    Index node;
    size_t edge;
};

/**
 * Implements two stacks which occupy the same space. This is only safe if
 * it is known that one is always smaller than the other, and that their
 * combined lengths never exceeds the capacity.
 *
 * Original author: David J. Pearce
 */
template<typename T>
struct DoubleStack {
    std::vector<T> items;
    size_t fp; // front pointer
    size_t bp; // back pointer

    DoubleStack(size_t capacity) :
        items(capacity), fp(0), bp(capacity)
    {}

    // ============================
    // Front stack
    // ============================

    bool empty_front() const {
        return fp == 0;
    }

    T& top_front() {
        return items[fp-1];
    }

    void pop_front() {
        fp--;
    }

    void push_front(const T& item) {
        items[fp++] = item;
    }

    // ============================
    // Back stack
    // ============================

    bool empty_back() const {
        return bp == items.size();
    }

    T& top_back() {
        return items[bp];
    }

    void pop_back() {
        bp++;
    }

    void push_back(const T& item) {
        items[--bp] = item;
	}
};

struct Transducer {
    size_t num_nodes;
    std::vector<std::vector<Edge>> succ_edges;

    /*
     * Create empty transducer
     */
    Transducer() :
        num_nodes(0)
    {}

    /*
     * Create identity transducer for given number of colors
     */
    Transducer(const Color colors) :
        num_nodes(1)
    {
        succ_edges.push_back({});
        for (Color c = 0; c < colors; c++) {
            add_edge(0, 0, c, c);
        }
    }

    /*
     * Create transducer for given number of color and set of tiles.
     */
    Transducer(const Color num_colors, const std::vector<Tile>& tiles) :
        num_nodes(num_colors)
    {
        succ_edges.resize(num_nodes);
        TileIndex t = 0;
        for (auto it = std::cbegin(tiles); it < std::cend(tiles); it++, t++) {
            const Tile& tile = *it;
            Edge* edge = add_edge(tile.west, tile.east, tile.north, tile.south);
            if (edge != nullptr) {
                edge->tiles.push_back(t);
            }
        }
    }

    Index add_node() {
        Index i = num_nodes;
        succ_edges.push_back({});
        num_nodes++;
        return i;
    }

    Edge* add_edge(Index source, Index target, Color north, Color south) {
        Edge edge(source, target, north, south);
        // do not add duplicate edges
        std::vector<Edge>& source_edges = succ_edges[source];
        const auto find_it = std::find(std::begin(source_edges), std::end(source_edges), edge);
        if (find_it == std::end(source_edges)) {
            succ_edges[source].push_back(std::move(edge));
            return &succ_edges[source].back();
        }
        else {
            return nullptr;
        }
    }

    inline void pearce_begin_visiting(
            std::vector<Index>& rindex,
            std::vector<bool>& root,
            DoubleStack<StackFrame>& stack,
            Index& index,
            Index v
    ) {
        stack.push_front({v, 0});
        root[v] = true;
        rindex[v] = index;
        index++;
    }

    void pearce_visit_loop(
            std::vector<Index>& rindex,
            std::vector<bool>& root,
            DoubleStack<StackFrame>& stack,
            Index& index,
            Index& c
    ) {
        StackFrame& frame = stack.top_front();
        const Index v = frame.node;
        const size_t end = succ_edges[v].size();
        for (size_t i = frame.edge; i <= end; i++) {
            if (i > 0) {
                // finish edge
                const Index w = succ_edges[v][i-1].target;
                if (rindex[w] < rindex[v]) {
                    rindex[v] = rindex[w];
                    root[v] = false;
                }
            }
            if (i < end) {
                // begin edge
                const Index w = succ_edges[v][i].target;
                if (rindex[w] == 0) {
                    frame.edge = i+1;
                    pearce_begin_visiting(rindex, root, stack, index, w);
                    return;
                }
            }
        }
        // finish visiting
        stack.pop_front();
        if (root[v]) {
            // create SCC
            index--;
            while (!stack.empty_back() && rindex[v] <= rindex[stack.top_back().node]) {
                Index w = stack.top_back().node;
                stack.pop_back();
                rindex[w] = c;
                index--;
            }
            rindex[v] = c;
            c--;
        }
        else {
            stack.push_back({v, 0});
        }
    }

    /* Pearce's algorithm for computing SCCs from
     * "A Space-Efficient Algorithm for Finding Strongly Connected Components"
     */
    std::vector<Index> compute_sccs() {
        std::vector<Index> rindex(num_nodes, 0);
        std::vector<bool> root(num_nodes, false);
        DoubleStack<StackFrame> stack(num_nodes);

        Index index = 1;
        Index c = num_nodes - 1;

        for (size_t v = 0; v < num_nodes; v++) {
            if (rindex[v] == 0) {
                pearce_begin_visiting(rindex, root, stack, index, v);
                while (!stack.empty_front()) {
                    pearce_visit_loop(rindex, root, stack, index, c);
                }
            }
        }

        return rindex;
    }

    /*
     * Simplify transducer by only keeping edges within the same strongly
     * connected component.
     * Note that we still keep nodes without edges here. Those will
     * however never be used in the next composition step.
     */
    void simplify() {
        std::vector<Index> scc_ids = compute_sccs();

        for (size_t i = 0; i < num_nodes; i++) {
            std::vector<Edge> new_succ_edges;
            for (size_t j = 0; j < succ_edges[i].size(); j++) {
                Edge& edge = succ_edges[i][j];
                if (scc_ids[edge.source] == scc_ids[edge.target]) {
                    new_succ_edges.push_back(std::move(edge));
                }
            }
            succ_edges[i].swap(new_succ_edges);
        }
    }

    inline void cycle_begin_visiting(
            std::vector<bool>& visited,
            std::vector<bool>& on_stack,
            std::stack<StackFrame, std::vector<StackFrame>>& stack,
            Index v
    ) {
        stack.push({v, 0});
        visited[v] = true;
        on_stack[v] = true;
    }

    bool cycle_visit_loop(
            std::vector<bool>& visited,
            std::vector<bool>& on_stack,
            std::stack<StackFrame, std::vector<StackFrame>>& stack,
            std::vector<Edge>& cycle
    ) {
        StackFrame& frame = stack.top();
        Index v = frame.node;
        const size_t end = succ_edges[v].size();
        for (size_t i = frame.edge; i < end; i++) {
            // begin edge
            const Edge& edge = succ_edges[v][i];
            if (edge.north == edge.south) {
                const Index w = edge.target;
                if (!visited[w]) {
                    // explore edge
                    frame.edge = i+1;
                    cycle_begin_visiting(visited, on_stack, stack, w);
                    return false;
                }
                else if (on_stack[w]) {
                    // found cycle
                    do {
                        cycle.push_back(succ_edges[v][i]);
                        stack.pop();
                        frame = stack.top();
                        v = frame.node;
                        i = frame.edge-1;
                    }
                    while (v != w);

                    std::reverse(std::begin(cycle), std::end(cycle));
                    return true;
                }
            }
        }
        // finish visiting
        stack.pop();
        on_stack[v] = false;
        return false;
    }

    std::vector<Edge> find_cycle() {
        std::vector<bool> visited(num_nodes, false);
        std::vector<bool> on_stack(num_nodes, false);

        std::vector<StackFrame> stack_container; stack_container.reserve(num_nodes);
        std::stack<StackFrame, std::vector<StackFrame>> stack(std::move(stack_container));

        std::vector<Edge> cycle;

        for (size_t v = 0; v < num_nodes; v++) {
            if (!visited[v]) {
                cycle_begin_visiting(visited, on_stack, stack, v);
                while (!stack.empty()) {
                    if (cycle_visit_loop(visited, on_stack, stack, cycle)) {
                        return cycle;
                    }
                }
            }
        }

        return cycle;
    }

    /*
     * Test if transducer contains a periodic cycle, i.e.
     * a cycle with same north and south color along each edge.
     * If the transducer contains a periodic cycle, it
     * returns the cycle as a rectangular tiling.
     */
    bool periodic(int height, Tiling& tiling) {
        // find cycle with dfs
        std::vector<Edge> cycle = find_cycle();

        if (!cycle.empty()) {
            int width = cycle.size();

            tiling.set_dimensions(width, height);
            for (int x = 0; x < width; x++) {
                std::vector<TileIndex>& column = cycle[x].tiles;
                for (int y = 0; y < height; y++) {
                    TileIndex i = column[y];
                    tiling.set_tile_index(x, y, i);
                }
            }
            return true;
        }
        return false;
    }

    bool empty() {
        for (size_t i = 0; i < num_nodes; i++) {
            if (!succ_edges[i].empty()) {
                return false;
            }
        }
        return true;
    }

    /*
     * Compose this transducer with another transducer,
     * and return the result.
     */
    Transducer compose(const Transducer& trans2) const {
        size_t product_size = num_nodes * trans2.num_nodes;
        // detect overflow
        assert (product_size / num_nodes == trans2.num_nodes);

        Transducer composition;
        std::vector<Index> node_map(product_size, -1);

        for (size_t i1 = 0; i1 < num_nodes; i1++) {
            for (const Edge& edge1 : succ_edges[i1]) {
                for (size_t i2 = 0; i2 < trans2.num_nodes; i2++) {
                    const Index i12 = i1 + i2*num_nodes;
                    for (const Edge& edge2 : trans2.succ_edges[i2]) {
                        if (edge1.south == edge2.north) {
                            Index j1 = edge1.target;
                            Index j2 = edge2.target;
                            const Index j12 = j1 + j2*num_nodes;

                            // add/find source node
                            Index pred_index = node_map[i12];
                            if (pred_index == -1) {
                                pred_index = composition.add_node();
                                node_map[i12] = pred_index;
                            }

                            // add/find target node
                            Index succ_index = node_map[j12];
                            if (succ_index == -1) {
                                succ_index = composition.add_node();
                                node_map[j12] = succ_index;
                            }

                            Edge* edge = composition.add_edge(pred_index, succ_index, edge1.north, edge2.south);
                            if (edge != nullptr) {
                                edge->tiles.reserve(edge1.tiles.size() + edge2.tiles.size());
                                edge->tiles.insert(std::end(edge->tiles), std::cbegin(edge1.tiles), std::cend(edge1.tiles));
                                edge->tiles.insert(std::end(edge->tiles), std::cbegin(edge2.tiles), std::cend(edge2.tiles));
                            }
                        }
                    }
                }
            }
        }

        return composition;
    }

    void print() const {
        for (size_t i = 0; i < num_nodes; i++) {
            for (const Edge& edge : succ_edges[i]) {
                std::cout << "  " << edge.source << " -> " << edge.target << " (" << (int)edge.north << "/" << (int)edge.south << ")" << std::endl;
            }
        }
    }
};

enum class TilesetClass { FINITE, PERIODIC, APERIODIC, UNKNOWN };

struct TilesetResult {
    TilesetClass result;
    Tiling tiling;
};

/*
 * Semi-decision algorithm to test if set of tiles is periodic or finite.
 * If it is periodic, also return a periodic tiling.
 */
TilesetResult test(const Tileset& tileset, int max_k, int verbosity = 0) {
    const Color num_colors = tileset.max_color + 1;
    Transducer trans(num_colors, tileset.tiles);
    if (verbosity >= 2) {
        std::cout << "Transducer for tileset without simplification" << std::endl;
        trans.print();
    }
    trans.simplify();
    Transducer trans_k(num_colors);

    TilesetResult r;
    Tiling tiling;
    for (int k = 1; k <= max_k; k++) {
        if (verbosity >= 1) {
            std::cout << "Testing k = " << k << std::endl;
        }
        trans_k = trans_k.compose(trans);
        if (verbosity >= 2) {
            std::cout << "Transducer for k = " << k << " before simplification" << std::endl;
            trans_k.print();
        }
        trans_k.simplify();
        if (verbosity >= 2) {
            std::cout << "Transducer for k = " << k << " after simplification" << std::endl;
            trans_k.print();
        }
        if (trans_k.empty()) {
            // tileset is finite
            r.result = TilesetClass::FINITE;
            return r;
        }
        if (trans_k.periodic(k, tiling)) {
            // tileset is periodic
            r.result = TilesetClass::PERIODIC;
            r.tiling = tiling;
            return r;
        }
    }
    // unknown result
    r.result = TilesetClass::UNKNOWN;
    return r;
}
