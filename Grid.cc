#include "Grid.h"
using Eigen::ArrayXd;

Grid::Grid(const size_t dim, const double lo, const double hi, 
    const size_t num_gc) : dimension(dim + 2 * num_gc), num_ghostcells(num_gc)
{
    dx = (hi - lo) / dim;
    low = lo - dx * num_ghostcells;
    high = hi + dx * num_ghostcells;
    fx = ArrayXd::LinSpaced(dimension + 1, low, high);
    cx = fx.segment(0, dimension) + dx / 2.0;
}

std::ostream &operator<< (std::ostream &os, const Grid &g)
{
    os << "Geometric limits: [" << g.low << ", " << g.high << "]" << std::endl;
    os << "Grid resolution: " << g.dimension << std::endl;
    os << "Grid space: dx = " << g.dx << std::endl;
    os << "Number of ghostcells: " << g.num_ghostcells << std::endl;
    return os;
}