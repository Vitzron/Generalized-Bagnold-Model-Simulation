#include "P_BC.h"

void P_BC::bc(Eigen::ArrayXd &phi, const size_t &num_gc) const
{
    Eigen::Index size = phi.size();

    for (size_t i = 0; i != num_gc; ++i)
    {
        phi(i) = phi(num_gc);
        phi(size - i - 1) = phi(size - num_gc - 1);
    }
}