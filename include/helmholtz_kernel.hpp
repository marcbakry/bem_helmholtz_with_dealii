#ifndef HELMHOLTZ_KERNEL_HPP
#define HELMHOLTZ_KERNEL_HPP

#include <cmath>

#include "defs_dealii.hpp"


class HelmholtzKernel
{
public:
    HelmholtzKernel(const std::complex<double> _kappa = std::complex<double>(0,0)): m_kappa(_kappa){};
    std::complex<double> single_layer(const dealii::Point<2> &_R);
    std::complex<double> single_layer(const dealii::Point<3> &_R);

    std::complex<double> fast_H1_0(const std::complex<double> _X);
    std::complex<double> fast_H1_0(const double _X);

private:

    std::complex<double> H1_0_low(const  std::complex<double> _X);
    std::complex<double> H1_0_high(const  std::complex<double> _X);

private:
    std::complex<double> m_kappa;
};

#endif