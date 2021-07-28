/*
 * File: Gas_EOS.h
 * Desc: header file of class Gas_EOS
 *       a isentropic gamma law for an ideal gas
 */

#ifndef Gas_EOS_
#define Gas_EOS_

#include <cmath>

class Gas_EOS
{
    public:
        Gas_EOS             () = default;
        ~Gas_EOS            () = default;
        Gas_EOS             (const double &p_0, const double &rho_0,
                             const double &gamma_0) : gamma(gamma_0), C(p_0 / std::pow(rho_0, gamma)) {};
        double p2den        (const double &p) const { return std::pow(p / C, 1.0 / gamma); };
    private:
        const double        gamma;
        const double        C;
};

#endif