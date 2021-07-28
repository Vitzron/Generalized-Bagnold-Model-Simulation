/*
 * File: U_BC.h
 * Desc: header file of class U_BC
 *       boundary condition on the velocity
 */

#ifndef U_BC_
#define U_BC_

#include <Eigen/Dense>

class U_BC
{
    public:
        U_BC        () = default;
        ~U_BC       () = default;
        void bc     (Eigen::ArrayXd &phi, const size_t &num_gc = 3) const;
};

#endif