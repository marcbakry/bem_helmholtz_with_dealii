#include "helmholtz_kernel.hpp"

std::complex<double> HelmholtzKernel::single_layer(const dealii::Point<2> &_R)
{
    auto kr = m_kappa*_R.norm();
    // compute using series
    return std::complex<double>(0,0.25)*fast_H1_0(kr);
}

std::complex<double> HelmholtzKernel::single_layer(const dealii::Point<3> &_R)
{
    auto r  = _R.norm();
    auto ikr = std::complex<double>(0,1)*m_kappa*r;
    return std::exp(ikr)/(r*4*dealii::numbers::PI);
}

std::complex<double> HelmholtzKernel::fast_H1_0(const std::complex<double> _X)
{ 
    return std::abs(1.0/3*_X) <= 1.0 ? H1_0_low(_X) : H1_0_high(_X);
}

std::complex<double> HelmholtzKernel::fast_H1_0(const double _X)
{
    return fast_H1_0(std::complex<double>(_X,0));
}

std::complex<double> HelmholtzKernel::H1_0_low(const std::complex<double> _X)
{
    auto T = _X*_X;
    auto J0 = std::complex<double>(0.00021,0);
    J0 = J0*T - 0.00029333;
    J0 = J0*T + 0.0444479;
    J0 = J0*T - 0.3163866;
    J0 = J0*T + 1.2656208;
    J0 = J0*T - 2.2499997;
    J0 = J0*T + 1.0;

    auto Y0 = -0.00024846*T + 0.00427916;
    Y0 = Y0*T - 0.04261214;
    Y0 = Y0*T + 0.25300117;
    Y0 = Y0*T - 0.74350384;
    Y0 = Y0*T + 0.60559366;
    Y0 = Y0*T + 0.36746691;
    Y0 += 2/dealii::numbers::PI*J0*std::log(0.5*_X);

    return J0 + std::complex<double>(0,1)*Y0;
}

std::complex<double> HelmholtzKernel::H1_0_high(const std::complex<double> _X)
{
    auto T = std::complex<double>(3,0)/_X;

    auto TH0 = 0.00013558*T - 0.00029333;
    TH0 = TH0*T - 0.00054125;
    TH0 = TH0*T + 0.00262573;
    TH0 = TH0*T - 0.00003954;
    TH0 = TH0*T - 0.04166397;
    TH0 = TH0*T - 0.78539816 + _X;

    auto F0 = 0.000144476*T - 0.00072805;
    F0 = F0*T + 0.00137237;
    F0 = F0*T - 0.00009512;
    F0 = F0*T - 0.00552740;
    F0 = F0*T - 0.00000077;
    F0 = F0*T + 0.79788456;

    return F0*std::exp(TH0)/std::sqrt(_X);
}