#include "GVF_BC.h"

void GVF_BC::bc(Eigen::ArrayXd &alpha) const
{
    Eigen::Index size = alpha.size();
    alpha(0) = alpha(1);
    alpha(size - 1) = alpha(size - 2);
}