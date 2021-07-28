/*
 * File: Eigen_Utilities.h
 * Desc: header file of class Eigen_Utilities
 *       added functions for Eigen::ArrayXXd
 */

#ifndef EIGEN_UTILITIES_
#define EIGEN_UTILITIES_

#include <Eigen/Dense>

class Eigen_Utilities
{
    public:
        Eigen_Utilities             () = default;
        ~Eigen_Utilities            () = default;
        void sign                   (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const;
        void positive               (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const;
        void negative               (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const
                                    { positive(oa, -ia); };
        void max                    (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia1, const Eigen::ArrayXXd &ia2) const;
        void max_zero_x             (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const;
        void min_zero_x             (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const;
        void minmod3                (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia1,
                                     const Eigen::ArrayXXd &ia2, const Eigen::ArrayXXd &ia3) const;
        void flipud                 (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const;
        void fliplr                 (Eigen::ArrayXXd &oa, const Eigen::ArrayXXd &ia) const;
        void vect2C_array           (Eigen::ArrayXXd &oa, const Eigen::VectorXd &v, 
                                     const Eigen::Index &rows, const Eigen::Index &cols) const;
        void vect2R_array           (Eigen::ArrayXXd &oa, const Eigen::VectorXd &v, 
                                     const Eigen::Index &rows, const Eigen::Index &cols) const;
        Eigen::ArrayXXd sign        (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd positive    (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd negative    (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd max         (const Eigen::ArrayXXd &ia1, const Eigen::ArrayXXd &ia2) const;
        Eigen::ArrayXXd max_zero_x  (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd min_zero_x  (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd flipud      (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd fliplr      (const Eigen::ArrayXXd &ia) const;
        Eigen::ArrayXXd minmod3     (const Eigen::ArrayXXd &ia1, const Eigen::ArrayXXd &ia2, const Eigen::ArrayXXd &ia3) const;
        Eigen::ArrayXXd vect2R_array(const Eigen::VectorXd &v, const Eigen::Index &rows, const Eigen::Index &cols) const;
        Eigen::ArrayXXd vect2C_array(const Eigen::VectorXd &v, const Eigen::Index &rows, const Eigen::Index &cols) const;
        Eigen::VectorXd array2R_vect(const Eigen::ArrayXXd &a) const;
        Eigen::VectorXd array2C_vect(const Eigen::ArrayXXd &a) const;
        Eigen::ArrayXXd central_difference_2rd_x(const Eigen::ArrayXXd &a, const double dx) const;
        Eigen::ArrayXXd central_difference_2rd_y(const Eigen::ArrayXXd &a, const double dy) const;
        Eigen::ArrayXXd second_partial_derivative_2rd_x(const Eigen::ArrayXXd &a, const double dx) const;
        Eigen::ArrayXXd second_partial_derivative_2rd_y(const Eigen::ArrayXXd &a, const double dy) const;
    private:
        const double                ZERO = 0.0;
        double sign                 (const double v) const;
	    double minmod3              (const double v1, const double v2, const double v3) const;
};

#endif