/*
 * File:    Mixture.h
 * Desc:    header file of class Mixture
 *          property of mixture of incompressible liquad and ideal gas
 */

#ifndef MIXTURE_
#define MIXTURE_

#include <Eigen/Dense>
#include "Gas_EOS.h"

class Mixture
{
    public:
        Mixture                 () = default;
        Mixture                 (const Gas_EOS &g, const double &rho) : a(g), rho_l(rho) {};
        ~Mixture                () = default;
        void mixture_density    (Eigen::ArrayXd &rho, const Eigen::ArrayXd &p,
                                 const Eigen::ArrayXd &alpha) const;
    private:
        Gas_EOS                 a;
        const double            rho_l;
};

#endif