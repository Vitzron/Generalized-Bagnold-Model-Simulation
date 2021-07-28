#include <iostream>
#include <vector>
#include "Pressure.h"
using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void Pressure::projection(Eigen::ArrayXd &p, const Eigen::ArrayXd &alpha,
        const Eigen::ArrayXd &rho, const Eigen::ArrayXd &velo,
        const double &dt, const double &dx) const
{
    if (p.size() != alpha.size() || p.size() != rho.size() || p.size() - 1 != velo.size())
    {
        cerr << "INCOMPATIBLE shapes in Pressure!" << endl;
        exit(EXIT_FAILURE);
    }

    Eigen::Index size = p.size() - 2;
    // right hand side
    ArrayXd sa;
    rhs(sa, alpha.segment(1, size), velo, dt, dx);
    // left hand side
    ArrayXd rho_x = rho.segment(2, size - 1) * rho.segment(1, size - 1) * (rho.segment(2, size - 1) + rho.segment(1, size - 1)) / (rho.segment(2, size - 1).square() + rho.segment(1, size - 1).square());
    SpMat A = SpMat(size, size);
    lhs(A, p.segment(1, size), alpha.segment(1, size), rho_x, dt, dx);
    // solver
    SimplicialLDLT<SpMat> solver;
    solver.compute(A);
    p.segment(1, size) = solver.solve(sa.matrix()).array();
}

void Pressure::rhs(Eigen::ArrayXd &sa, 
        const Eigen::ArrayXd &alpha, const Eigen::ArrayXd &velo,
        const double &dt, const double &dx) const
{
    Eigen::Index size = alpha.size();
    sa = (velo.segment(1, size) - velo.segment(0, size)) / dx / dt - alpha / gamma / dt / dt;
}

void Pressure::lhs(Eigen::SparseMatrix<double> &A,
        const Eigen::ArrayXd &p, const Eigen::ArrayXd &alpha,
        const Eigen::ArrayXd &rho,
        const double &dt, const double &dx) const
{
    Eigen::Index size = p.size();
    std::vector<T> tripletList;
    tripletList.reserve(size * 3);
    double X = 1.0 / dx / dx;

    // first row
    tripletList.push_back(T(0, 0, -alpha(0) / gamma / p(0) / dt / dt - X / rho(0)));
    tripletList.push_back(T(0, 1, X / rho(0)));

    // middle rows
    for (Index i = 1; i != size - 1; ++i)
    {
        tripletList.push_back(T(i, i - 1, X / rho(i - 1)));
        tripletList.push_back(T(i, i + 1, X / rho(i)));
        tripletList.push_back(T(i, i, -alpha(i) / gamma / p(i) / dt / dt - X / rho(i - 1) - X / rho(i)));
    }

    // last row
    tripletList.push_back(T(size - 1, size - 2, X / rho(size - 2)));
    tripletList.push_back(T(size - 1, size - 1, -alpha(size - 1) / gamma / p(size - 1) / dt / dt - X / rho(size - 2)));

    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A.makeCompressed();
}