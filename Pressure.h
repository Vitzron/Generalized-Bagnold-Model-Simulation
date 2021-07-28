/*
 * File: Pressure.h
 * Desc: header file of class Pressure
 */

#ifndef PRESSURE_
#define PRESSURE_

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Pressure
{
    public:
        Pressure        () = default;
        Pressure        (const double &gamma_0) : gamma(gamma_0) {};
        ~Pressure       () = default;
        void projection (Eigen::ArrayXd &p, const Eigen::ArrayXd &alpha,
                         const Eigen::ArrayXd &rho, const Eigen::ArrayXd &velo,
                         const double &dt, const double &dx) const;
    private:
        void lhs        (Eigen::SparseMatrix<double> &A,
                         const Eigen::ArrayXd &p, const Eigen::ArrayXd &alpha,
                         const Eigen::ArrayXd &rho,
                         const double &dt, const double &dx) const;
        void rhs        (Eigen::ArrayXd &sa, 
                         const Eigen::ArrayXd &alpha, const Eigen::ArrayXd &velo,
                         const double &dt, const double &dx) const;
        const double    gamma;
};

#endif