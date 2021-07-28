#include "THINC.h"
#include <iostream>
using namespace std;
using namespace Eigen;

constexpr double beta = 3.5;
constexpr double eps = 1.0e-10;

void THINC::flux(double &f, const Eigen::Array3d &phi, 
        const double &u, const double &dt, const double &dx) const
{
    double gamma, lambda;
    Index iup;
    double x_m;

    if (u == 0.0)
    {
        f = 0.0;
    }
    else
    {
        iup = u < 0.0 ? 2 : 1;
        
        if (fabs(phi(iup)) < eps || (fabs(phi(iup) - 1.0)) < eps)
        {
            f = phi(iup) * u * dt;
        }
        else
        {
            gamma = phi(0) < phi(2) ? 1.0 : -1.0;
            lambda = u < 0.0 ? 0.0 : 1.0;
            x_m = 0.5 / beta * log((exp(beta / gamma * (1.0 + gamma - 2.0 * phi(iup))) - 1.0) 
                / (1.0 - exp(beta / gamma * (1.0 - gamma - 2.0 * phi(iup)))));
            f = 0.5 * (u * dt - gamma * dx / beta * log((cosh(beta * (lambda - x_m - u * dt / dx))) 
                / (cosh(beta * (lambda - x_m)))));
        }
    }
}

void THINC::update(Eigen::ArrayXd &phi, const Eigen::ArrayXd &velo,
        const double &dt, const double &dx) const
{
    if (phi.size() + 1 != velo.size())
    {
        cerr << "INCOMPATIBLE sizes!" << endl;
        exit(EXIT_FAILURE);
    }
    if (phi.maxCoeff()> 1.0 + eps || phi.minCoeff() < -eps)
    {
        cerr << "WRONG values!" << endl;
        exit(EXIT_FAILURE);
    }

    Index size = phi.size() - 1;
    ArrayXd phi_flux = ArrayXd::Zero(size);
    phi_flux(0) = 0.0;
    phi_flux(size - 1) = 0.0;
    for (Index i = 0; i != size - 2; ++i)
    {
        flux(phi_flux(i + 1), phi.segment(i, 3), velo(i + 2), dt, dx);
    }
    phi.segment(1, size - 1) = phi.segment(1, size - 1) 
        - (phi_flux.segment(1, size - 1) - phi_flux.segment(0, size - 1)) / dx
        + phi.segment(1, size - 1) * (velo.segment(2, size - 1) - velo.segment(1, size - 1)) * dt / dx;
}