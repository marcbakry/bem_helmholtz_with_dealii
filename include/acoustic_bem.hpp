#ifndef ACOUSTIC_BEM_HPP
#define ACOUSTIC_BEM_HPP

#include <fstream>

#include "defs_dealii.hpp"
#include "helmholtz_kernel.hpp"

class AcousticBEM
{
public:
    AcousticBEM(const std::complex<double> _kappa);

    void run();

private:
    void setup_grids();
    void setup_system();
    void assemble_matrices();
    void solve();
    void write_outputs() const;

    bool isCloseInteraction(const dealii::Triangulation<1,2>::cell_iterator &_cell1,const dealii::Triangulation<1,2>::cell_iterator &_cell2) const;

private:
    bool m_write_mesh;
    unsigned int m_qorder_singular;
    unsigned int m_qorder_close;
    unsigned int m_qorder_far;
    double m_close_coef; // coefficient permettant de determiner si une interaction est proche ou lointaine

    // physique
    HelmholtzKernel m_kernel;
    std::complex<double> m_kappa;

    // BEM
    dealii::Triangulation<1,2> m_triangulation;
    dealii::FE_Q<1,2> m_fe;
    dealii::DoFHandler<1,2> m_dof_handler;
    dealii::MappingQ<1,2> m_mapping;
    dealii::FullMatrix<std::complex<double>> m_lhs;
    dealii::Vector<std::complex<double>> m_rhs;
    dealii::Vector<std::complex<double>> m_sol;
};

#endif