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
#include <vector>
#include <string>

void print_usage(const char* argv[]) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] tileset [output.svg]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "    -h          Show this help" << std::endl;
    std::cout << "    -k <max_k>  Maximum value of k to construct transducer for" << std::endl;
    std::cout << "    -r          Randomize colors in output" << std::endl;
    std::cout << "    -a          Always print largest horizontal tiling found, even if non-periodic vertically" << std::endl;
    std::cout << "    -t          Only test for periodic/finite tilings, do not keep information to recover tiling" << std::endl;
    std::cout << "    -v          Increase verbosity (may be specified more than once)" << std::endl;

}

void run_test(
        const Tileset& tileset,
        const int max_k,
        const bool print_always,
        const bool create_tiling,
        const bool randomize_colors,
        const bool output_svg,
        const std::string& svg_filename,
        const int verbosity
) {
    TilesetResult r = test(tileset, max_k, print_always, create_tiling, verbosity);

    if (r.result == TilesetClass::FINITE) {
        std::cout << "finite" << std::endl;
    }
    else if (r.result == TilesetClass::PERIODIC) {
        std::cout << "periodic" << std::endl;
    }
    else if (r.result == TilesetClass::APERIODIC) {
        // this will currently never happen
        std::cout << "aperiodic" << std::endl;
    }
    else {
        std::cout << "unknown" << std::endl;
    }

    if (create_tiling && (r.result == TilesetClass::PERIODIC || print_always)) {
        r.tiling.print();
        if (output_svg) {
            draw_tiling(tileset, r.tiling, svg_filename, randomize_colors);
        }
    }
}

int main(const int argc, const char* argv[]) {
    if (argc < 2) {
        print_usage(argv);
        return EXIT_SUCCESS;
    }

    // parse arguments
    std::string tileset_filename;
    std::string svg_filename;
    bool output_svg = false;
    int max_k = 30;
    int verbosity = 0;

    bool parse_input_name = true;
    bool parse_output_name = true;
    bool parse_k = false;
    bool randomize_colors = false;
    bool print_always = false;
    bool create_tiling = true;

    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (!arg.compare(0, 1, "-")) {
            if (!arg.compare(1, 1, "k")) {
                parse_k = true;
            }
            else if (!arg.compare(1, 1, "v")) {
                verbosity++;
            }
            else if (!arg.compare(1, 1, "h")) {
                print_usage(argv);
                return EXIT_SUCCESS;
            }
            else if (!arg.compare(1, 1, "r")) {
                randomize_colors = true;
            }
            else if (!arg.compare(1, 1, "a")) {
                print_always = true;
            }
            else if (!arg.compare(1, 1, "t")) {
                create_tiling = false;
            }
            else {
                std::cerr << "Error: unknown option: " << arg << std::endl;
                return EXIT_FAILURE;
            }
        }
        else if (parse_k) {
            parse_k = false;
            try {
                max_k = std::stoi(arg);
            }
            catch(std::exception const & e) {
                std::cerr << "Error : could not parse argument for '-k': " << arg << std::endl;
                return EXIT_FAILURE;
            }
        }
        else if (parse_input_name) {
            parse_input_name = false;
            tileset_filename = arg;
        }
        else if (parse_output_name) {
            parse_output_name = false;
            svg_filename = arg;
            output_svg = true;
        }
        else {
            std::cerr << "Error: additional argument: " << arg << std::endl;
            return EXIT_FAILURE;
        }
    }
    if (tileset_filename.empty()) {
        std::cerr << "Error: no input file given" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream tileset_file;
    tileset_file.open(tileset_filename, std::ios_base::in);
    if (tileset_file.fail()) {
        std::cerr << "Error: Failed to open " << tileset_filename << std::endl;
        return EXIT_FAILURE;
    }

    try {
        const Tileset tileset = Tileset::parse_tileset(tileset_file);
        run_test(tileset, max_k, print_always, create_tiling, randomize_colors, output_svg, svg_filename, verbosity);
        return EXIT_SUCCESS;
    }
    catch(const std::bad_alloc&) {
        std::cerr << "Error: out of memory" << std::endl;
        return EXIT_FAILURE;
    }
    catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "Error: unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
