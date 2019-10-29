/*
 * Implemantation of the semi-decision algorithm described in to
 * test if a set of tiles is aperiodic, using transducers.
 *
 * Described in:
 * Emmanuel Jeandel and Michael Rao: An aperiodic set of 11 Wang tiles
 *
 * (optimizations for minimization and bisimulation not implemented)
 */

#include "tiles.hpp"
#include "transducer.hpp"
#include "visualize.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

int main(const int argc, const char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " tileset [output.svg]" << std::endl;
        return EXIT_SUCCESS;
    }
    std::string tileset_filename(argv[1]);
    std::string svg_filename;
    bool output_svg = false;
    if (argc >= 3) {
        svg_filename = std::string(argv[2]);
        output_svg = true;
    }

    std::ifstream tileset;
    tileset.open(tileset_filename, std::ios_base::in);
    if (tileset.fail()) {
        std::cerr << "Error: Failed to open " << tileset_filename << std::endl;
        return EXIT_FAILURE;
    }

    constexpr int MAX_K = 30;

    std::vector<Tile> tiles;
    int max_color = 0;

    std::string line;
    while (std::getline(tileset, line)) {
        // check if line is not a comment
        if (line.compare(0, 1, "#")) {
            std::stringstream linestream(line);
            int north, south, west, east;
            linestream >> north >> south >> west >> east;
            tiles.push_back(Tile(north, south, west, east));
            max_color = std::max(max_color, north);
            max_color = std::max(max_color, south);
            max_color = std::max(max_color, west);
            max_color = std::max(max_color, east);
        }
    }
    const int num_colors = max_color + 1;

    Tiling tiling;
    int result = test_periodic(tiles, num_colors, MAX_K, tiling);
    if (result == 0) {
        std::cout << "finite" << std::endl;
    }
    else if (result == 1) {
        std::cout << "periodic " << tiling.get_width() << " " << tiling.get_height() << std::endl;
        tiling.print();
        if (output_svg) {
            draw_tiling(tiles, tiling, svg_filename);
        }
    }
    else {
        std::cout << "unknown" << std::endl;
    }

    return EXIT_SUCCESS;
}
