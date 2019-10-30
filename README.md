# Tiling

Implementation of a semi-decision procedure for the [tiling/domino problem](https://en.wikipedia.org/wiki/Wang_tile) based on transducers as described in the paper "[An aperiodic set of 11 Wang tiles](https://arxiv.org/abs/1506.06492)" by Emmanuel Jeandel and Michael Rao.

## Problem statement

The problem is to decide whether a set of [Wang tiles](https://en.wikipedia.org/wiki/Wang_tile), i.e. square tiles defined by their colors on the north/south/west/east side, can tile the infinite plane.
Moreover, one can state the problem as classifying a set of tiles as one the following:

- Periodic: the tileset can tile the plane with a periodic pattern, i.e. they can be composed into a finite repeatable rectangular pattern.
- Aperiodic: the tileset can tile the plane, but not with a periodic pattern.
- Finite: the tileset can not tile the plane, and can therefore also not tile a finite square of some size.

The problem is undecidable, however deciding if a set of Wang tiles is periodic or finite is semi-decidable.

## Building

Simply invoke `make` to build the program. You need GCC in a version supporting the C++14 standard, e.g. GCC 8.0.

## Usage

The program can be used as follows:
```
bin/tiling [OPTIONS] tileset [output.svg]
```
Here, `tileset` is the file name of the tileset to be analyzed. A few examples are given in the `tilesets` folder.
With `output.svg` one can specify an output file for a periodic pattern, if one is found.
The options include `-k MAX_k` to specify the maximum value of `k` until which a transducer in the composition is constructed.
Invoke the program with the `-h` option for a list of other options.

### Input

An input file describes a tileset by listing on each line a tile in the format `n s w e` with non-negative integers `n`, `s`, `w`, `e`, describing the colors on the north, south, west and east side of the tile, respectively.
The input can also contain line comments starting with `#`.

### Output

The output is either the word `finite`, if the tileset is deduced to be finite, `periodic` if the tileset is deduced to be periodic, or `unknown` if the maximum bound for composing the transducer is reached.
In the periodic case, additionally a rectangular periodic pattern is printed  giving in each cell the index of the tile used.
The pattern can also be printed in SVG format by additionally giving an output file as an argument.

## References

- Emmanuel Jeandel and Michael Rao: An aperiodic set of 11
Wang tiles
  - [Paper on arXiv](https://arxiv.org/abs/1506.06492)
  - [Implementation by the authors](https://framagit.org/mrao/small-wang-tile-sets)
- Emmanuel Jeandel and Pascal Vanier: The Undecidability of the Domino Problem
  - [Lecture notes](https://www.lacl.fr/~pvanier/rech/cirm.pdf)
  - Slides: [Part 1](https://www.cirm-math.fr/ProgWeebly/Renc1720/Jeandel.pdf),
    [Part 2](https://www.cirm-math.fr/ProgWeebly/Renc1720/Jeandel2.pdf),
    [Part 3](https://www.cirm-math.fr/ProgWeebly/Renc1720/Jeandel3.pdf)
