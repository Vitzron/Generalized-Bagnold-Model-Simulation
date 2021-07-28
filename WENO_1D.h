/*
 * File: WENO_1D.h
 * Desc: header file of class WENO_1D
 *       one-dimensional weighted ENO
 */

#ifndef WENO_1D_
#define WENO_1D_

#include <Eigen/Dense>

class WENO_1D
{
    public:
        WENO_1D         () = default;
        ~WENO_1D        () = default;
        void update     (Eigen::ArrayXd &phi, const Eigen::ArrayXd &velo,
                         const double &dt, const double &dx, const size_t &num_gc = 3) const;
    private:
        void set_phi_s  (Eigen::ArrayXd &phi_s, const Eigen::ArrayXd &v1,
                         const Eigen::ArrayXd &v2, const Eigen::ArrayXd &v3,
                         const Eigen::ArrayXd &v4, const Eigen::ArrayXd &v5) const;
};

#endif