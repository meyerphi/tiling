#pragma once

#include <vector>

#include "tiles.hpp"
#include "simple_svg.hpp"

// colors selected from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
constexpr int colorset[][3] = {
    {230, 25, 75}, // red
    {255, 225, 25}, // yellow
    {60, 180, 75}, // green
    {0, 130, 200}, // blue
    {240, 50, 230}, // magenta
    {70, 240, 240}, // cyan
    {245, 130, 48}, // orange
    {230, 190, 255}, // lavender
    {170, 110, 40}, // brown
    {0, 128, 128}, // teal
    {128, 0, 0}, // maroon
    {0, 0, 128}, // navy
    {250, 190, 190}, // pink
    {255, 250, 200}, // beige
    {170, 255, 195}, // mint
    {145, 30, 180}, // purple
    {255, 215, 180}, // apricot
    {128, 128, 0}, // olive
    {210, 245, 60}, // lime
    {128, 128, 128}, // grey
    {0, 0, 0}, // black
    {255, 255, 255} // white
};
template <typename T, size_t N>
constexpr size_t dim(T(&)[N]) { return N; }
constexpr size_t num_available_colors = dim(colorset);

svg::Color get_color(const size_t i) {
    return svg::Color(colorset[i][0], colorset[i][1], colorset[i][2]);
}
svg::Stroke get_stroke() {
    return svg::Stroke(0.025, svg::Color(0,0,0));
}

void add_tile(const Tile& tile, svg::Document& doc, int x, int y) {
    svg::Polygon north(get_color(tile.north), get_stroke());
    svg::Polygon east(get_color(tile.east), get_stroke());
    svg::Polygon south(get_color(tile.south), get_stroke());
    svg::Polygon west(get_color(tile.west), get_stroke());

    north << svg::Point(x,y) << svg::Point(x+1,y) << svg::Point(x+0.5,y+0.5);
    east << svg::Point(x+1,y) << svg::Point(x+1,y+1) << svg::Point(x+0.5,y+0.5);
    south << svg::Point(x+1,y+1) << svg::Point(x,y+1) << svg::Point(x+0.5,y+0.5);
    west << svg::Point(x,y+1) << svg::Point(x,y) << svg::Point(x+0.5,y+0.5);

    svg::Rectangle border(svg::Point(x,y), 1, 1, svg::Fill(), get_stroke());

    doc << north << east << south << west << border;
}

bool draw_tiling(const Tileset& tileset, const Tiling& tiling, std::string filename) {
    if (tileset.max_color >= num_available_colors) {
        std::cerr << "Error: not enough colors available to draw tiling" << std::endl;
        return false;
    }

    svg::Dimensions dimensions(tiling.get_width(), tiling.get_height());
    svg::Document doc(filename, svg::Layout(dimensions, svg::Layout::TopLeft));

    for (int x = 0; x < tiling.get_width(); x++) {
        for (int y = 0; y < tiling.get_height(); y++) {
            const TileIndex i = tiling.get_tile_index(x, y);
            const Tile& tile = tileset.get_tile(i);

            add_tile(tile, doc, x, y);
        }
    }

    doc.save();
    return true;
}
