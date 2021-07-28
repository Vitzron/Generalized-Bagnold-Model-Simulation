/*
 * File: GVF_BC.h
 * Desc: header file of class GVF_BC
 *       boundary condition on the volume fraction of gas
 */

#ifndef GVF_BC_
#define GVF_BC_

#include <Eigen/Dense>

class GVF_BC
{
    public:
        GVF_BC      () = default;
        ~GVF_BC     () = default;
        void bc     (Eigen::ArrayXd &alpha) const;
};

#endif