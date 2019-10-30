#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>

typedef uint8_t TileIndex;
typedef uint8_t Color;

struct Tile {
    const Color north;
    const Color south;
    const Color west; // source in transducer
    const Color east; // target in transducer

    Tile(Color north, Color south, Color west, Color east) :
        north(north), south(south), west(west), east(east)
    {}

    bool operator==(const Tile& rhs) const {
        return north == rhs.north && south == rhs.south &&
               west  == rhs.west  && east  == rhs.east;
    }
    bool operator!=(const Tile& rhs) const {
        return !operator==(rhs);
    }
};

struct Tileset {
public:
    const Color max_color;
    const std::vector<Tile> tiles;

private:
    Tileset(Color mc, std::vector<Tile> t) :
        max_color(mc),
        tiles(std::move(t))
    {}

public:
    const Tile& get_tile(TileIndex i) const {
        assert(i >= 0);
        assert(i < tiles.size());
        return tiles[i];
    }

    static Tileset parse_tileset(std::istream& in) {
        int max_color = 0;
        std::vector<Tile> tiles;

        std::string line;
        while (std::getline(in, line)) {
            // check if line is not empty or a a comment
            if (!line.empty() && line.compare(0, 1, "#")) {
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

        return Tileset(max_color, std::move(tiles));
    }
};

struct Tiling {
private:
    int _width;
    int _height;
    std::vector<TileIndex> tile_indices;

public:
    Tiling() : _width(0), _height(0) {}

    int get_width() const {
        return _width;
    }

    int get_height() const {
        return _height;
    }

    void set_dimensions(int width, int height) {
        _width = width;
        _height = height;
        tile_indices.resize(width*height);
    }

    TileIndex get_tile_index(int x, int y) const {
        assert(x >= 0);
        assert(y >= 0);
        assert(x < _width);
        assert(y < _height);
        return tile_indices[x + _width*y];
    }

    void set_tile_index(int x, int y, TileIndex i) {
        assert(x >= 0);
        assert(y >= 0);
        assert(x < _width);
        assert(y < _height);
        tile_indices[x + _width*y] = i;
    }

    void print() const {
        for (int y = 0; y < _height; y++) {
            for (int x = 0; x < _width; x++) {
                std::cout << (int)get_tile_index(x, y);
                if (x + 1 < _width) {
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        }
    }
};
