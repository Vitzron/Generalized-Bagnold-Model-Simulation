#include <iostream>
#include <fstream>
#include "Grid.h"
#include "THINC.h"
#include "WENO_1D.h"
#include "Pressure.h"
#include "GVF_BC.h"
#include "U_BC.h"
#include "P_BC.h"
#include "Gas_EOS.h"
#include "Mixture.h"
#include "Eigen_Utilities.h"
using namespace std;
using namespace Eigen;

constexpr double g = 9.80;
constexpr double bar = 1.0e5;

int main(int argc, char* argv[])
{
    const double rho_l = atof(argv[1]);
    const double rho_g = atof(argv[2]);
    const double p0 = atof(argv[3]);
    const double gamma = atof(argv[4]);
    const double scale = atof(argv[5]);
    const char *bottom_pressure = argv[6];

    const size_t dim = 300;
    const size_t num_gc = 3;

    double lo = 0.0 / scale;
    double hi = 15.0 / scale;
    Grid grid(dim, lo, hi);
    double dx = grid.get_dx();
    cout << grid << endl;

    ArrayXd fx = grid.get_fx();
    ArrayXd cx = grid.get_cx();
    // cout << "fx: " << endl << fx.transpose() << endl;
    // cout << "cx: " << endl << cx.transpose() << endl;

    ArrayXd alpha = ArrayXd::Zero(dim + 2);
    ArrayXd u = ArrayXd::Zero(dim + 2 * num_gc - 1);
    ArrayXd p = ArrayXd::Zero(dim + 2 * num_gc);

    double lo_l = 2.0;
    double hi_l = 10.0;

    lo_l /= scale;
    hi_l /= scale;

    Mixture mix(Gas_EOS(p0, rho_g, gamma), rho_l);
    THINC thinc = THINC();
    WENO_1D weno = WENO_1D();
    Pressure pressure(gamma);
    U_BC u_bc = U_BC();
    P_BC p_bc = P_BC();
    GVF_BC g_bc = GVF_BC();

    Eigen_Utilities eu = Eigen_Utilities();
    alpha = 1.0 - eu.positive((cx.segment(2, dim + 2) - lo_l) * (hi_l - cx.segment(2, dim + 2)));
    // cout << "alpha: " << endl << alpha.transpose() << endl;
    p = p0;
    // cout << "p: " << endl << p.transpose() << endl;

    const double dt_0 = 0.0001;
    double dt = dt_0 / scale;
    cout << "dx: " << dx << "; dt: " << dt << endl;

    ofstream output;
    output.open(bottom_pressure, std::ios::out);
    if (!output)
    {
        cerr << "CANNOT open output file" << bottom_pressure << "!" << endl;
        exit(EXIT_FAILURE);
    }
    output << 0.0 << " " << p0 / bar << endl;

    size_t steps = 2000000;
    const double T_total = 1.3;
    const size_t records = 1300;
    const size_t count = ceil(ceil(T_total / sqrt(scale) / dt) / records);
    double T = 0.0;

    ArrayXd velo;
    ArrayXd u_t, p_t;
    ArrayXd alpha_0, u_0, p_0;
    ArrayXd rho, rho_x;
    ArrayXd div;

    for (size_t i = 0; i != steps; ++i)
    {
        cout << "step: " << i << "; T: " << T * sqrt(scale) << "; bottom pressure: " << (p(num_gc) - p0) * scale << "; volume of water: " << (alpha.segment(1, dim).size() - alpha.segment(1, dim).sum()) * dx;
        cout << "; max of alpha: " << alpha.maxCoeff() << "; min of alpha: " << alpha.minCoeff() << endl;

        alpha_0 = alpha;
        u_0 = u;
        p_0 = p;

        for (size_t substep = 0; substep != 3; ++substep)
        {
            u_t = u;

            //*********************advection*********************
            // air volume fraction
            velo = u_t.segment(1, dim + 3);
            thinc.update(alpha, velo, dt, dx);
            g_bc.bc(alpha);

            // u
            velo = u_t.segment(3, dim - 1);
            weno.update(u, velo, dt, dx);
            u_bc.bc(u);

            // p
            velo = (u_t.segment(2, dim) + u_t.segment(3, dim)) / 2.0;
            weno.update(p, velo, dt, dx);
            p_bc.bc(p);

            //*********************nonadvection (i)*********************
            u -= g * dt;
            u_bc.bc(u);

            //*********************nonadvection (ii)********************
            velo = u.segment(2, dim + 1);
            mix.mixture_density(rho, p.segment(2, dim + 2), alpha);
            rho_x = rho.segment(1, dim - 1) * rho.segment(2, dim - 1) * (rho.segment(1, dim - 1) + rho.segment(2, dim - 1)) / (rho.segment(1, dim - 1).square() + rho.segment(2, dim - 1).square());
            p_t = p.segment(2, dim + 2);
            pressure.projection(p_t, alpha, rho, velo, dt, dx);
            p.segment(2, dim + 2) = p_t;
            p_bc.bc(p);

            // update velocity
            u.segment(3, dim - 1) -= dt / rho_x * (p.segment(4, dim - 1) - p.segment(3, dim - 1)) / dx;
            u_bc.bc(u);

            velo = u.segment(2, dim + 1);
            div = (velo.segment(1, dim) - velo.segment(0, dim)) / dx;
            alpha.segment(1, dim) += (dt * div * (1.0 - alpha.segment(1, dim)));

            if (substep == 1)
            {
                alpha = 3.0 / 4.0 * alpha_0 + 1.0 / 4.0 * alpha;
                u = 3.0 / 4.0 * u_0 + 1.0 / 4.0 * u;
                p = 3.0 / 4.0 * p_0 + 1.0 / 4.0 * p;
            }
            if (substep == 2)
            {
                alpha = 1.0 / 3.0 * alpha_0 + 2.0 / 3.0 * alpha;
                u = 1.0 / 3.0 * u_0 + 2.0 / 3.0 * u;
                p = 1.0 / 3.0 * p_0 + 2.0 / 3.0 * p;
            }
        }

        velo = u.segment(2, dim + 1);
        div = (velo.segment(1, dim) - velo.segment(0, dim)) / dx;
        cout << "Error index: " << (div * (1.0 - alpha).segment(1, dim)).abs().sum() << endl;
        if (i % count == 0)
            output << T * sqrt(scale) << " " << (p(3) - p0) / bar * scale + 1.0 << endl;

        T += dt;
        if (T * sqrt(scale) > T_total)
            break;
    }

    output.flush();
    output.close();

    return 0;
}