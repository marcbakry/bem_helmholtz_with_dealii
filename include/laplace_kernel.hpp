#ifndef LAPLACE_KERNEL_HPP
#define LAPLACE_KERNEL_HPP

#include <cmath>

#include "defs_dealii.hpp"

class LaplaceKernel
{
public:
    LaplaceKernel(){};
    double single_layer(const dealii::Tensor<1,2> &_R);
    double single_layer(const dealii::Tensor<1,3> &_R);
};

#endif