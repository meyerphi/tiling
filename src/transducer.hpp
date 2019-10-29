#pragma once

#include <vector>
#include <stack>
#include <map>
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
};

struct Transducer {
    size_t num_nodes;
    std::vector<Edge> edges;
    std::vector<std::vector<Index>> succ_edges;

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
        succ_edges.resize(1);
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
            // do not add duplicate tiles
            const auto find_it = std::find(std::cbegin(tiles), it, tile);
            if (find_it == it) {
                Edge& edge = add_edge(tile.west, tile.east, tile.north, tile.south);
                // only keep index of tile in edge
                edge.tiles.push_back(t);
            }
        }
    }

    Index add_node() {
        Index i = num_nodes;
        succ_edges.push_back({});
        num_nodes++;
        return i;
    }

    Edge& add_edge(Index source, Index target, Color north, Color south) {
        Index e = edges.size();
        edges.push_back(Edge(source, target, north, south));
        succ_edges[source].push_back(e);
        return edges[e];
    }

    /*
     * Tarjan's algorithm for computing SCCs.
     */
    void compute_sccs(std::vector<Index>& scc_ids) {
        struct StackFrame {
            Index node;
            size_t edge;
        };

        std::stack<StackFrame> call_stack;
        std::vector<bool> on_stack(num_nodes, false);
        std::vector<Index> index(num_nodes, -1);
        std::vector<Index> lowlink(num_nodes, -1);
        std::stack<Index> stack;

        Index cur_scc_id = 0;

        size_t i = 0;
        Index cur_index = 0;
        bool new_node = false;
        bool done = false;

        Index v;
        size_t j;
        while (!done) {
            if (new_node) {
                index[v] = cur_index;
                lowlink[v] = cur_index;
                cur_index++;
                stack.push(v);
                on_stack[v] = true;
                j = 0;
                new_node = false;
            }
            else {
                if (call_stack.empty()) {
                    while (i < num_nodes && index[i] != -1) {
                        i++;
                    }
                    if (i == num_nodes) {
                        done = true;
                    }
                    v = i;
                    new_node = true;
                }
                else {
                    StackFrame frame = call_stack.top();
                    call_stack.pop();
                    v = frame.node;
                    j = frame.edge;
                    Index w = edges[succ_edges[v][j]].target;
                    lowlink[v] = std::min(lowlink[v], lowlink[w]);
                    j++;
                }
            }
            while (!new_node && j != succ_edges[v].size()) {
                Index e = succ_edges[v][j];
                Edge& edge = edges[e];
                Index w = edge.target;

                if (index[w] == -1) {
                    call_stack.push({v, j});
                    v = w;
                    new_node = true;
                }
                else if (on_stack[w]) {
                    lowlink[v] = std::min(lowlink[v], index[w]);
                }
                j++;
            }
            if (!new_node) {
                if (lowlink[v] == index[v]) {
                    Index w;
                    do {
                        w = stack.top();
                        stack.pop();
                        scc_ids[w] = cur_scc_id;
                        on_stack[w] = false;
                    }
                    while (v != w);
                    cur_scc_id++;
                }
            }
        }
    }

    /*
     * Simplify transducer by only keeping edges within the same strongly
     * connected component.
     * Note that we still keep nodes without edges here. Those will
     * however never be used in the next composition step.
     */
    void simplify() {
        std::vector<Index> scc_ids(num_nodes, -1);
        compute_sccs(scc_ids);

        for (size_t i = 0; i < num_nodes; i++) {
            std::vector<Index> new_succ_edges;
            for (size_t j = 0; j < succ_edges[i].size(); j++) {
                Index e = succ_edges[i][j];
                Edge& edge = edges[e];
                Index v = edge.source;
                Index w = edge.target;
                if (scc_ids[v] == scc_ids[w]) {
                    new_succ_edges.push_back(e);
                }
            }
            succ_edges[i].swap(new_succ_edges);
        }
    }

    /*
     * DFS to find cycle with same north and south
     * colors along each edge.
     */
    bool find_cycle(std::vector<Index>& cycle) {
        struct StackFrame {
            Index node;
            size_t edge;
        };

        std::stack<StackFrame> stack;
        std::vector<Index> parent(num_nodes, -1);
        std::vector<bool> on_stack(num_nodes, false);
        std::vector<Index> index(num_nodes, -1);

        size_t i = 0;
        Index cycle_node = -1;
        Index cur_index = 0;
        bool new_node = false;
        bool done = false;

        Index v = -1;
        size_t j = 0;
        while (!done) {
            if (new_node) {
                // explore new node
                on_stack[v] = true;
                index[v] = cur_index;
                cur_index++;
                j = 0;
                new_node = false;
            }
            else {
                if (stack.empty()) {
                    // find next initial node
                    while (i < num_nodes && index[i] != -1) {
                        i++;
                    }
                    if (i == num_nodes) {
                        done = true;
                    }
                    else {
                        v = i;
                        new_node = true;
                    }
                }
                else {
                    // explore next edge after backtracking
                    StackFrame frame = stack.top();
                    stack.pop();
                    v = frame.node;
                    j = frame.edge;
                    j++;
                }
            }
            while (!new_node && j != succ_edges[v].size()) {
                Index e = succ_edges[v][j];
                Edge& edge = edges[e];
                if (edge.north == edge.south) {
                    Index w = edge.target;
                    if (index[w] == -1) {
                        // forward edge, explore
                        parent[w] = e;
                        stack.push({v, j});
                        v = w;
                        new_node = true;
                    }
                    else if (on_stack[w]) {
                        // backedge, found cycle
                        parent[w] = e;
                        cycle_node = w;
                        new_node = true;
                        done = true;
                    }
                }
                j++;
            }
            if (!new_node) {
                // backtrack from current node
                on_stack[v] = false;
            }
        }

        if (cycle_node >= 0) {
            Index v = cycle_node;
            do {
                Index e = parent[v];
                v = edges[e].source;
                cycle.push_back(e);
            }
            while (v != cycle_node);
            std::reverse(std::begin(cycle), std::end(cycle));
            return true;
        }
        else {
            return false;
        }
    }

    /*
     * Test if transducer contains a periodic cycle, i.e.
     * a cycle with same north and south color along each edge.
     * If the transducer contains a periodic cycle, it
     * returns the cycle as a rectangular tiling.
     */
    bool periodic(int height, Tiling& tiling) {
        // find cycle with dfs
        std::vector<Index> cycle;

        if (find_cycle(cycle)) {
            int width = cycle.size();

            tiling.set_dimensions(width, height);
            for (int x = 0; x < width; x++) {
                std::vector<TileIndex>& column = edges[cycle[x]].tiles;
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
        Transducer composition;
        std::map<std::pair<Index, Index>, Index> node_map;

        for (size_t i1 = 0; i1 < num_nodes; i1++) {
            for (Index e1 : succ_edges[i1]) {
                for (size_t i2 = 0; i2 < trans2.num_nodes; i2++) {
                    for (Index e2 : trans2.succ_edges[i2]) {
                        const Edge& edge1 = edges[e1];
                        const Edge& edge2 = trans2.edges[e2];

                        if (edge1.south == edge2.north) {
                            Index j1 = edge1.target;
                            Index j2 = edge2.target;

                            // add/find source node
                            auto pred_it = node_map.find({i1,i2});
                            Index pred_index;
                            if (pred_it == std::end(node_map)) {
                                pred_index = composition.add_node();
                                node_map.insert({{i1, i2}, pred_index});
                            }
                            else {
                                pred_index = pred_it->second;
                            }

                            // add/find target node
                            auto succ_it = node_map.find({j1,j2});
                            Index succ_index;
                            if (succ_it == std::end(node_map)) {
                                succ_index = composition.add_node();
                                node_map.insert({{j1, j2}, succ_index});
                            }
                            else {
                                succ_index = succ_it->second;
                            }

                            Edge& e = composition.add_edge(pred_index, succ_index, edge1.north, edge2.south);
                            e.tiles.reserve(edge1.tiles.size() + edge2.tiles.size());
                            e.tiles.insert(std::end(e.tiles), std::cbegin(edge1.tiles), std::cend(edge1.tiles));
                            e.tiles.insert(std::end(e.tiles), std::cbegin(edge2.tiles), std::cend(edge2.tiles));
                        }
                    }
                }
            }
        }

        return composition;
    }

    void print() const {
        std::cout << "Transducer with n = " << num_nodes << "; m = " << edges.size() << std::endl;
        for (size_t i = 0; i < num_nodes; i++) {
            std::cout << "Node " << i << ":";
            std::cout << std::endl;
            for (Index e : succ_edges[i]) {
                std::cout << "  " << e << " = " << edges[e].source << " -> " << edges[e].target << "(" << (int)edges[e].north << "/" << (int)edges[e].south << ")" << std::endl;
            }
            std::cout << std::endl;
        }
    }
};

/*
 * Semi-decision algorithm to test if set of tiles is periodic or finite.
 */
int test_periodic(const std::vector<Tile>& tiles, Color num_colors, int max_k, Tiling& tiling) {
    Transducer trans(num_colors, tiles);
    Transducer trans_k(num_colors);
    trans.simplify();

    for (int k = 1; k <= max_k; k++) {
        trans_k = trans_k.compose(trans);
        trans_k.simplify();
        //std::cout << "Transducer for k = " << k << " after simplification" << std::endl;
        //trans_k.print();
        if (trans_k.empty()) {
            // tileset is finite
            return 0;
        }
        if (trans_k.periodic(k, tiling)) {
            // tileset is periodic
            return 1;
        }
    }
    // unknown result
    return -1;
}
