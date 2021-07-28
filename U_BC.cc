#include "U_BC.h"

void U_BC::bc(Eigen::ArrayXd &phi, const size_t &num_gc) const
{
    Eigen::Index size = phi.size();

    for (size_t i = 0; i != num_gc - 1; ++i)
    {
        phi(i) = -phi(2 * num_gc - i - 2);
        phi(size - i - 1) = -phi(size - 2 * num_gc + i + 1);
    }
    phi(num_gc - 1) = 0.0;
    phi(size - num_gc) = 0.0;
}