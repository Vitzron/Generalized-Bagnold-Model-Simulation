#include <iostream>
#include "Mixture.h"
using namespace std;
using namespace Eigen;

void Mixture::mixture_density(Eigen::ArrayXd &rho, const Eigen::ArrayXd &p, 
        const Eigen::ArrayXd &alpha) const
{
    if (p.size() != alpha.size())
    {
        cerr << "INCOMPATIBLE sizes in Mixture!" << endl;
        exit(EXIT_FAILURE);
    }

    Eigen::Index size = p.size();
    rho = ArrayXd::Zero(size);

    for (Eigen::Index i = 0; i != size; ++i)
    {
        rho(i) = a.p2den(p(i)) * alpha(i) + rho_l * (1.0 - alpha(i));
    }
}