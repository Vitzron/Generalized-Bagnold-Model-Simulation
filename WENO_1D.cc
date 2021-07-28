#include <iostream>
#include "WENO_1D.h"
#include "Eigen_Utilities.h"
using namespace std;
using namespace Eigen;

constexpr double SMALL = 1.0e-99;
constexpr double small = 1.0e-6;

void WENO_1D::update(Eigen::ArrayXd &phi, const Eigen::ArrayXd &velo,
        const double &dt, const double &dx, const size_t &num_gc) const
{
    Eigen::Index size = phi.size() - 2 * num_gc;

    if (size != velo.size())
    {
        cerr << "INCOMPATIBLE shapes in WENO_1D!" << endl;
        exit(EXIT_FAILURE);
    }

    Eigen_Utilities eu = Eigen_Utilities();
    ArrayXd v, phi_x_m, phi_x_p;

    // phi_x_minus
    v = (phi.segment(1, phi.size() - 2) - phi.segment(0, phi.size() - 2)) / dx;
    set_phi_s(phi_x_m, v.segment(0, size), v.segment(1, size), v.segment(2, size), v.segment(3, size), v.segment(4, size));
    // phi_x_plus
    v = (phi.segment(2, phi.size() - 2) - phi.segment(1, phi.size() - 2)) / dx;
    set_phi_s(phi_x_p, v.segment(4, size), v.segment(3, size), v.segment(2, size), v.segment(1, size), v.segment(0, size));

    phi.segment(num_gc, size) -= dt * (eu.max_zero_x(velo) * phi_x_m + eu.min_zero_x(velo) * phi_x_p);
}

void WENO_1D::set_phi_s(Eigen::ArrayXd &phi_s, const Eigen::ArrayXd &v1,
        const Eigen::ArrayXd &v2, const Eigen::ArrayXd &v3,
        const Eigen::ArrayXd &v4, const Eigen::ArrayXd &v5) const
{
    ArrayXd phi_x_1 = v1 / 3.0 - v2 * 7.0 / 6.0 + v3 * 11.0 / 6.0;
    ArrayXd phi_x_2 = -v2 / 6.0 + v3 * 5.0 / 6.0 + v4 / 3.0;
    ArrayXd phi_x_3 = v3 / 3.0 + v4 * 5.0 / 6.0 - v5 / 6.0;

    ArrayXd S1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3).square() 
        + 1.0 / 4.0 * (v1 - 4.0 * v2 + 3.0 * v3).square();
    ArrayXd S2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4).square()
        + 1.0 / 4.0 * (v2 - v4).square();
    ArrayXd S3 = 13.0 / 12.0 * (v3 - 2.0 * v4 + v5).square()
        + 1.0 / 4.0 * (3.0 * v3 - 4.0 * v4 + v5).square();
    
    Eigen_Utilities eu = Eigen_Utilities();
    ArrayXd epsilon = small * eu.max(eu.max(eu.max(eu.max(v1.square(), v2.square()), v3.square()), v4.square()), v5.square()) + SMALL;

    ArrayXd alpha1 = 0.1 / (S1 + epsilon).square();
    ArrayXd alpha2 = 0.6 / (S2 + epsilon).square();
    ArrayXd alpha3 = 0.3 / (S3 + epsilon).square();
    ArrayXd alpha = alpha1 + alpha2 + alpha3;
    phi_s = alpha1 / alpha * phi_x_1 + alpha2 / alpha * phi_x_2 + alpha3 / alpha * phi_x_3;
}