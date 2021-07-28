/*
 * File:    THINC.h
 * Desc:    header file of classs THINC
 *          angent of hyperbola for interface capturing
 */

#ifndef THINC_
#define THINC_

#include <Eigen/Dense>

class THINC
{
    public:
        THINC               () = default;
        ~THINC              () = default;
        void update         (Eigen::ArrayXd &phi, const Eigen::ArrayXd &velo,
                             const double &dt, const double &dx) const;
    private:
        void flux           (double &f, const Eigen::Array3d &phi,
                             const double &u, const double &dt, const double &dx) const;
};

#endif