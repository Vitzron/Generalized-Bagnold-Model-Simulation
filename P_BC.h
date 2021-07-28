/*
 * File: P_BC.h
 * Desc: header file of class P_BC
 */

#ifndef P_BC_
#define P_BC_

#include <Eigen/Dense>

class P_BC
{
    public:
        P_BC        () = default;
        ~P_BC       () = default;
        void bc     (Eigen::ArrayXd &phi, const size_t &num_gc = 3) const;
};

#endif