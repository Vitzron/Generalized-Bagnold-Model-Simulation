/*
 * File: Grid.h
 * Desc: header file of class Grid
 */

#ifndef GRID_
#define GRID_

#include <iostream>
#include <Eigen/Dense>

class Grid
{
    friend std::ostream &operator<< (std::ostream &os, const Grid &g);
    public:
        Grid                    () = default;
        Grid                    (const size_t dim, const double lo, 
                                 const double hi, const size_t num_gc = 3);
        ~Grid                   () = default;
        size_t get_numgc        () const { return num_ghostcells; };
        size_t get_dim          () const { return dimension; };
        double get_lo           () const { return low; };
        double get_hi           () const { return high; };
        double get_dx           () const { return dx; };
        Eigen::ArrayXd get_fx   () const { return fx; };
        Eigen::ArrayXd get_cx   () const { return cx; };
    private:
        size_t                  dimension;
        size_t                  num_ghostcells;
        double                  low;
        double                  high;
        double                  dx;
        Eigen::ArrayXd          fx;
        Eigen::ArrayXd          cx;
};

#endif