#include <iostream>
#include "Eigen_Utilities.h"
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using std::cerr;
using std::endl;

// ************************************************************************************************
// public functions
// ************************************************************************************************
void Eigen_Utilities::sign(ArrayXXd &oa, const ArrayXXd &ia) const
{
    oa = ia;
    for (Eigen::Index j = 0; j != ia.cols(); ++j)
    {
        for (Eigen::Index i = 0; i != ia.rows(); ++i)
        {
            oa(i, j) = sign(ia(i, j));
        }
    }
}

void Eigen_Utilities::positive(ArrayXXd &oa, const ArrayXXd &ia) const
{
    oa = ia;
    for (Eigen::Index j = 0; j != ia.cols(); ++j)
    {
        for (Eigen::Index i = 0; i != ia.rows(); ++i)
        {
            oa(i, j) = ia(i, j) > ZERO ? 1.0 : ZERO;
        }
    }
}

void Eigen_Utilities::max(ArrayXXd &oa, const ArrayXXd &ia1, const ArrayXXd &ia2) const
{
    if (ia1.rows() != ia2.rows() || ia1.cols() != ia2.cols())
    {
        cerr << "INCOMPATIBLE shapes!" << endl;
        exit(EXIT_FAILURE);
    }
    oa = ia1;
    for (Eigen::Index j = 0; j != ia1.cols(); ++j)
    {
        for (Eigen::Index i = 0; i != ia1.rows(); ++i)
        {
            oa(i, j) = fmax(ia1(i, j), ia2(i, j));
        }
    }
}

void Eigen_Utilities::max_zero_x(ArrayXXd &oa, const ArrayXXd &ia) const
{
    ArrayXXd pa;
    positive(pa, ia);
    oa = pa * ia;
}

void Eigen_Utilities::min_zero_x(ArrayXXd &oa, const ArrayXXd &ia) const
{
    ArrayXXd na;
    negative(na, ia);
    oa = na * ia;
}

void Eigen_Utilities::minmod3(ArrayXXd &oa,
        const ArrayXXd &ia1, const ArrayXXd &ia2, const ArrayXXd &ia3) const
{
    if (ia1.rows() != ia2.rows() ||
        ia1.cols() != ia2.cols() ||
        ia2.rows() != ia3.rows() ||
        ia2.cols() != ia3.cols())
    {
        cerr << "INCOMPATIBLE shapes!" << endl;
        exit(EXIT_FAILURE);
    }
    oa = ia1;
    for (Eigen::Index j = 0; j != ia1.cols(); ++j)
    {
        for (Eigen::Index i = 0; i != ia1.rows(); ++i)
        {
            oa(i, j) = minmod3(ia1(i, j), ia2(i, j), ia3(i, j));
        }
    }
}

void Eigen_Utilities::flipud(ArrayXXd &oa, const ArrayXXd &ia) const
{
    oa = ia;
    for (Eigen::Index i = 0; i != ia.rows(); ++i)
    {
        oa.row(i) = ia.row(ia.rows() - 1 - i);
    }
}

void Eigen_Utilities::fliplr(ArrayXXd &oa, const ArrayXXd &ia) const
{
    oa = ia;
    for (Eigen::Index j = 0; j != ia.cols(); ++j)
    {
        oa.col(j) = ia.col(ia.cols() - 1 - j);
    }
}

void Eigen_Utilities::vect2C_array(ArrayXXd &oa,
        const VectorXd &v, const Eigen::Index &rows, const Eigen::Index &cols) const
{
    if (v.size() < rows * cols)
    {
        cerr << "INCOMPATIBLE shapes!" << endl;
        exit(EXIT_FAILURE);
    }
    oa = ArrayXXd::Zero(rows, cols);
    for (Eigen::Index j = 0; j != cols; ++j)
    {
        oa.col(j) = v.segment(j * rows, rows);
    }
}

void Eigen_Utilities::vect2R_array(ArrayXXd &oa,
        const VectorXd &v, const Eigen::Index &rows, const Eigen::Index &cols) const
{
    if (v.size() < rows * cols)
    {
        cerr << "INCOMPATIBLE shapes!" << endl;
        exit(EXIT_FAILURE);
    }
    oa = ArrayXXd::Zero(rows, cols);
    for (Eigen::Index i = 0; i != rows; ++i)
    {
        oa.row(i) = v.segment(i * cols, cols);
    }
}

ArrayXXd Eigen_Utilities::sign(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    sign(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::positive(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    positive(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::negative(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    negative(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::max(const ArrayXXd &ia1, const ArrayXXd &ia2) const
{
    ArrayXXd oa = ia1;
    max(oa, ia1, ia2);
    return oa;
}

ArrayXXd Eigen_Utilities::max_zero_x(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    max_zero_x(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::min_zero_x(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    min_zero_x(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::fliplr(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    fliplr(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::flipud(const ArrayXXd &ia) const
{
    ArrayXXd oa = ia;
    flipud(oa, ia);
    return oa;
}

ArrayXXd Eigen_Utilities::minmod3(const ArrayXXd &ia1, const ArrayXXd &ia2, const ArrayXXd &ia3) const
{
    ArrayXXd oa = ia1;
    minmod3(oa, ia1, ia2, ia3);
    return oa;
}

ArrayXXd Eigen_Utilities::vect2R_array(const VectorXd &v, const Eigen::Index &rows, const Eigen::Index &cols) const
{
    ArrayXXd oa(rows, cols);
    vect2R_array(oa, v, rows, cols);
    return oa;
}

ArrayXXd Eigen_Utilities::vect2C_array(const VectorXd &v, const Eigen::Index &rows, const Eigen::Index &cols) const
{
    ArrayXXd oa(rows, cols);
    vect2C_array(oa, v, rows, cols);
    return oa;
}

VectorXd Eigen_Utilities::array2C_vect(const ArrayXXd &a) const
{
    VectorXd v(a.size());
    for (Eigen::Index j = 0; j != a.cols(); ++j)
    {
        v.segment(j * a.rows(), a.rows()) = a.col(j);
    }
    return v;
}

VectorXd Eigen_Utilities::array2R_vect(const ArrayXXd &a) const
{
    VectorXd v(a.size());
    for (Eigen::Index i = 0; i != a.rows(); ++i)
    {
        v.segment(i * a.cols(), a.cols()) = a.row(i);
    }
    return v;
}

ArrayXXd Eigen_Utilities::central_difference_2rd_x(const ArrayXXd &a, const double dx) const
{
    return (a.rightCols(a.cols() - 2) - a.leftCols(a.cols() - 2)) / 2.0 / dx;
}

ArrayXXd Eigen_Utilities::central_difference_2rd_y(const ArrayXXd &a, const double dy) const
{
    return (a.topRows(a.rows() - 2) - a.bottomRows(a.rows() - 2)) / 2.0 / dy;
}

ArrayXXd Eigen_Utilities::second_partial_derivative_2rd_x(const ArrayXXd &a, const double dx) const
{
    return (a.rightCols(a.cols() - 2) + a.leftCols(a.cols() - 2)
        - 2.0 * a.middleCols(1, a.cols() - 2)) / dx / dx;
}

ArrayXXd Eigen_Utilities::second_partial_derivative_2rd_y(const ArrayXXd &a, const double dy) const
{
    return (a.topRows(a.rows() - 2) + a.bottomRows(a.rows() - 2)
        - 2.0 * a.middleRows(1, a.rows() - 2)) / dy / dy;
}

// ************************************************************************************************
// private functions
// ************************************************************************************************
double Eigen_Utilities::sign(const double v) const
{
    if (v == ZERO)
        return ZERO;
    else
        return v > ZERO ? 1.0 : -1.0;
}

double Eigen_Utilities::minmod3(const double v1, const double v2, const double v3) const
{
    if (sign(v1) != sign(v2) ||
        sign(v1) != sign(v3) ||
        sign(v2) != sign(v3))
        return ZERO;
    else
        return sign(v1) * fmin(fmin(fabs(v1), fabs(v2)), fabs(v3));
}