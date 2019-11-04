#pragma once

#include <vector>
#include <stack>
#include <limits>
#include <algorithm>

#include "tiles.hpp"

typedef unsigned int Index;
constexpr Index UNEXPLORED = std::numeric_limits<Index>::max();

struct Edge {
    Index target;
    Color north;
    Color south;
    std::vector<TileIndex> tiles;

    // only needed to resize edge array to smaller array
    // will never be called
    Edge () { assert(false); }

    Edge(Index target, Color north, Color south) :
        target(target), north(north), south(south)
    {}

    /*
     * Note: equality function does not compare the set of tiles,
     * as it is sufficient to keep one edge with otherwise same values.
     */
    bool operator==(const Edge& rhs) const {
        return target == rhs.target && north == rhs.north && south == rhs.south;
    }
    bool operator!=(const Edge& rhs) const {
        return !operator==(rhs);
    }
};

// optimized stack frame for Tarjan and DFS
struct StackFrame {
    Index node;
    Index edge;
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
    std::vector<std::vector<Edge>> succ_edges;
    std::vector<Index> nodes;
    std::vector<Edge> edges;

    /*
     * Create empty transducer
     */
    Transducer() {
        nodes.push_back(0);
    }

    /*
     * Create identity transducer for given number of colors
     */
    Transducer(const Color colors) : Transducer() {
        Index i = add_node();
        for (Color c = 0; c < colors; c++) {
            add_edge(i, i, c, c);
        }
    }

    /*
     * Create transducer for given number of color and set of tiles.
     */
    Transducer(const Color num_colors, const std::vector<Tile>& tiles) : Transducer() {
        // group tiles by source color
        std::vector<std::vector<std::pair<Tile, TileIndex>>> source_tiles(num_colors);

        TileIndex t = 0;
        for (const Tile& tile : tiles) {
            source_tiles[tile.west].push_back({tile, t});
            t++;
        }

        for (const auto& tileset : source_tiles) {
            Index i = add_node();
            for (const auto& entry : tileset) {
                const Tile& tile = entry.first;
                const TileIndex t = entry.second;
                assert(i == tile.west);
                Edge& edge = add_edge(tile.west, tile.east, tile.north, tile.south);
                edge.tiles.push_back(t);
            }
        }
    }

    Index num_nodes() const {
        return nodes.size()-1;
    }

    Index num_edges() const {
        return edges.size();
    }

    Index get_edge_begin(Index v) const {
        return nodes[v];
    }
    Index get_edge_end(Index v) const {
        return nodes[v+1];
    }
    const Edge& get_edge(Index e) const {
        return edges[e];
    }

    Index add_node() {
        Index i = num_nodes();
        nodes[i] = num_edges();
        nodes.push_back(num_edges());
        return i;
    }

    Edge& add_edge(Index source, Index target, Color north, Color south) {
        assert(nodes[source+1] == edges.size());

        nodes[source+1]++;
        edges.emplace_back(target, north, south);
        return edges.back();
    }

    inline void pearce_begin_visiting(
            std::vector<Index>& rindex,
            std::vector<bool>& root,
            DoubleStack<StackFrame>& stack,
            Index& index,
            Index v
    ) const {
        stack.push_front({v, get_edge_begin(v)});
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
    ) const {
        StackFrame& frame = stack.top_front();
        const Index v = frame.node;
        const Index begin = get_edge_begin(v);
        const Index end = get_edge_end(v);
        for (Index i = frame.edge; i <= end; i++) {
            if (i > begin) {
                // finish edge
                const Index w = get_edge(i-1).target;
                if (rindex[w] < rindex[v]) {
                    rindex[v] = rindex[w];
                    root[v] = false;
                }
            }
            if (i < end) {
                // begin edge
                const Index w = get_edge(i).target;
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
    std::vector<Index> compute_sccs() const {
        std::vector<Index> rindex(num_nodes(), 0);
        std::vector<bool> root(num_nodes(), false);
        DoubleStack<StackFrame> stack(num_nodes());

        Index index = 1;
        Index c = num_nodes() - 1;

        for (Index v = 0; v < num_nodes(); v++) {
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

        Index cur_start = 0;
        Index cur_index = 0;

        for (Index i = 0; i < num_nodes(); i++) {
            const Index source_scc_id = scc_ids[i];
            for (Index j = get_edge_begin(i); j < get_edge_end(i); j++) {
                Edge& edge = edges[j];
                if (source_scc_id == scc_ids[edge.target]) {
                    // apparently copy is here faster than move
                    // with guarded check cur_index
                    edges[cur_index] = edge;
                    cur_index++;
                }
            }
            nodes[i] = cur_start;
            cur_start = cur_index;
        }
        nodes[num_nodes()] = cur_start;
        edges.resize(cur_start);
    }

    inline void cycle_begin_visiting(
            std::vector<bool>& visited,
            std::vector<bool>& on_stack,
            std::stack<StackFrame, std::vector<StackFrame>>& stack,
            Index v
    ) const {
        stack.push({v, get_edge_begin(v)});
        visited[v] = true;
        on_stack[v] = true;
    }

    bool cycle_visit_loop(
            std::vector<bool>& visited,
            std::vector<bool>& on_stack,
            std::stack<StackFrame, std::vector<StackFrame>>& stack,
            std::vector<Edge>& cycle,
            const bool periodic
    ) const {
        StackFrame& frame = stack.top();
        const Index v = frame.node;
        const Index end = get_edge_end(v);
        for (Index i = frame.edge; i < end; i++) {
            // begin edge
            const Edge& edge = get_edge(i);
            if (!periodic || edge.north == edge.south) {
                const Index w = edge.target;
                if (!visited[w]) {
                    // explore edge
                    frame.edge = i+1;
                    cycle_begin_visiting(visited, on_stack, stack, w);
                    return false;
                }
                else if (on_stack[w]) {
                    // found cycle
                    cycle.push_back(get_edge(i));
                    while (stack.top().node != w) {
                        stack.pop();
                        frame = stack.top();
                        cycle.push_back(get_edge(frame.edge-1));
                    }
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

    void print_debug() {
        std::cout << "Node array:" << std::endl;
        for (Index i = 0; i < num_nodes(); i++) {
            std::cout << " [" << nodes[i] << "," << nodes[i+1] << "]";
        }
        std::cout << std::endl;
        std::cout << "Edge array:" << std::endl;
        for (const Edge& e : edges) {
            std::cout << " " << e.target;
        }
        std::cout << std::endl;
    }

    std::vector<Edge> find_cycle(const bool periodic) const {

        std::vector<bool> visited(num_nodes(), false);
        std::vector<bool> on_stack(num_nodes(), false);

        std::vector<StackFrame> stack_container; stack_container.reserve(num_nodes());
        std::stack<StackFrame, std::vector<StackFrame>> stack(std::move(stack_container));

        std::vector<Edge> cycle;

        for (Index v = 0; v < num_nodes(); v++) {
            if (!visited[v]) {
                cycle_begin_visiting(visited, on_stack, stack, v);
                while (!stack.empty()) {
                    if (cycle_visit_loop(visited, on_stack, stack, cycle, periodic)) {
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
    bool periodic(int height, Tiling& tiling, const bool proper) const {
        // find cycle with dfs
        std::vector<Edge> cycle = find_cycle(proper);

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

    bool empty() const {
        for (Index i = 0; i < num_nodes(); i++) {
            if (get_edge_begin(i) != get_edge_end(i)) {
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
        const Index product_size = num_nodes() * trans2.num_nodes();
        // detect overflow
        assert (product_size / num_nodes() == trans2.num_nodes());

        Transducer composition;
        std::vector<Index> node_map(product_size, UNEXPLORED);

        for (Index i1 = 0; i1 < num_nodes(); i1++) {
            for (Index i2 = 0; i2 < trans2.num_nodes(); i2++) {
                const Index i12 = i1 + i2*num_nodes();
                Index source = UNEXPLORED;
                for (Index k1 = get_edge_begin(i1); k1 < get_edge_end(i1); k1++) {
                    const Edge& edge1 = get_edge(k1);
                    for (Index k2 = trans2.get_edge_begin(i2); k2 < trans2.get_edge_end(i2); k2++) {
                        const Edge& edge2 = trans2.get_edge(k2);
                        if (edge1.south == edge2.north) {

                            const Index j1 = edge1.target;
                            const Index j2 = edge2.target;
                            const Index target = j1 + j2*num_nodes();
                            if ((j1 < i1 || (j1 == i1 && j2 < i2)) && node_map[target] == UNEXPLORED) {
                                // target node visited but not added, can skip edge
                                continue;
                            }

                            if (source == UNEXPLORED) {
                                source = composition.add_node();
                                node_map[i12] = source;
                            }

                            Edge& edge = composition.add_edge(source, target, edge1.north, edge2.south);
                            // this part can be omitted if one only needs to decide periodicity,
                            // but is needed to actually to recover a periodic tiling
                            edge.tiles.reserve(edge1.tiles.size() + edge2.tiles.size());
                            edge.tiles.insert(std::end(edge.tiles), std::cbegin(edge1.tiles), std::cend(edge1.tiles));
                            edge.tiles.insert(std::end(edge.tiles), std::cbegin(edge2.tiles), std::cend(edge2.tiles));
                        }
                    }
                }
            }
        }

        // set correct edge targets
        Index sink = composition.add_node();

        for (Index e = 0; e < composition.num_edges(); e++) {
            Edge& edge = composition.edges[e];
            Index map_target = node_map[edge.target];
            if (map_target == UNEXPLORED) {
                // target node not added, edge will be removed by simplify
                map_target = sink;
            }
            edge.target = map_target;
        }

        return composition;
    }

    void print_size() const {
        std::cout << "n = " << num_nodes() << "; m = " << num_edges() << std::endl;
    }

    void print() const {
        for (Index i = 0; i < num_nodes(); i++) {
            for (Index j = get_edge_begin(i); j < get_edge_end(i); j++) {
                const Edge& edge = get_edge(j);
                std::cout << "  " << i << " -> " << edge.target << " (" << (int)edge.north << "/" << (int)edge.south << ")";
                std::cout << " [";
                for (const TileIndex& t : edge.tiles) {
                    std::cout << " " << (int)t;
                }
                std::cout << " ]" << std::endl;
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
TilesetResult test(const Tileset& tileset, int max_k, bool always_test = false, int verbosity = 0) {
    const Color num_colors = tileset.max_color + 1;
    Transducer trans(num_colors, tileset.tiles);
    if (verbosity >= 2) {
        std::cout << "Transducer for tileset without simplification" << std::endl;
        trans.print();
        trans.print_debug();
    }
    trans.simplify();
    Transducer trans_k(num_colors);

    TilesetResult r;
    for (int k = 1; k <= max_k; k++) {
        if (verbosity >= 1) {
            std::cout << "Testing k = " << k << std::endl;
        }
        trans_k = trans_k.compose(trans);
        if (verbosity >= 1) {
            std::cout << "Size of composition: ";
            trans_k.print_size();
        }
        if (verbosity >= 2) {
            std::cout << "Transducer for k = " << k << " before simplification" << std::endl;
            trans_k.print();
            trans_k.print_debug();
        }
        trans_k.simplify();
        if (verbosity >= 2) {
            std::cout << "Transducer for k = " << k << " after simplification" << std::endl;
            trans_k.print();
            trans_k.print_debug();
        }
        if (trans_k.empty()) {
            // tileset is finite
            r.result = TilesetClass::FINITE;
            return r;
        }
        if (trans_k.periodic(k, r.tiling, true)) {
            // tileset is periodic
            r.result = TilesetClass::PERIODIC;
            return r;
        }
        if (always_test) {
            trans_k.periodic(k, r.tiling, false);
        }
    }
    // unknown result
    r.result = TilesetClass::UNKNOWN;
    return r;
}
