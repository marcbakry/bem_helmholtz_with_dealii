#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>

#include "defs_dealii.hpp"
#include "helmholtz_kernel.hpp"


void test_hankel_functions();

int main()
{
    //
    test_hankel_functions();
    //
    return EXIT_SUCCESS;
}

void test_hankel_functions()
{
    // initialize arguments
    unsigned int N = 100;
    double xMax = 10.0;
    std::vector<double> x(N+1);
    auto dX = xMax/N;
    for(auto i=0; i<x.size(); ++i) x[i] = i*dX;
    // bessel functions
    std::vector<double> J(x.size()), Y(x.size());
    std::vector<std::complex<double>> H(x.size());
    HelmholtzKernel m_kernel(1.0);

    for(auto i=0; i<J.size(); ++i)
    {
        //
        J[i] = std::cyl_bessel_j(0,x[i]);
        Y[i] = std::cyl_neumann(0,x[i]);
        H[i] = m_kernel.fast_H1_0(x[i]);
    }

    // output to file
    std::ofstream file("../data/comparaison_hankel.csv");
    for(auto i=0; i<x.size(); ++i)
    {
        file << x[i] << ";" << J[i] << ";" << Y[i] << ";" << H[i].real() << ";" << H[i].imag() << std::endl;
    }
    file.close();
}