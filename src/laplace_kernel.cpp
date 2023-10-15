#include "laplace_kernel.hpp"

double LaplaceKernel::single_layer(const dealii::Tensor<1,2> &_R)
{
    return -std::log(_R.norm())/(2*dealii::numbers::PI);
}

double LaplaceKernel::single_layer(const dealii::Tensor<1,3> &_R)
{
    return 1.0/(_R.norm()*4*dealii::numbers::PI);
}