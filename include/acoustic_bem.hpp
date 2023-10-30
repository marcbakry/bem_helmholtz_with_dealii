#ifndef ACOUSTIC_BEM_HPP
#define ACOUSTIC_BEM_HPP

#include <fstream>
#include <cassert>
#include <complex>
#include <cmath>
#include <utility>
#include <cassert>

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
    void radiate();
    void write_outputs() const;

    unsigned int isCloseInteraction(const dealii::Triangulation<1,2>::cell_iterator &_cell1,const dealii::Triangulation<1,2>::cell_iterator &_cell2) const;

    void buildNonSingularCellMatrix(const dealii::Triangulation<1,2>::cell_iterator &_tcell,const dealii::Triangulation<1,2>::cell_iterator &_scell,dealii::FEValues<1,2> &_tfe, dealii::FEValues<1,2> &_sfe, dealii::FullMatrix<std::complex<double>> &_cmatrix);

    void buildSingularCellMatrix(const dealii::Triangulation<1,2>::cell_iterator &_tcell,const dealii::Triangulation<1,2>::cell_iterator &_scell,dealii::FEValues<1,2> &_tfe,dealii::FullMatrix<std::complex<double>> &_cmatrix);

    void buildRHSCellVector(const dealii::Triangulation<1,2>::cell_iterator &_tcell,dealii::FEValues<1,2> &_tfe,dealii::Vector<std::complex<double>> &_cvector);

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
    dealii::LAPACKFullMatrix<std::complex<double>> m_lhs;
    dealii::Vector<std::complex<double>> m_rhs;
    dealii::Vector<std::complex<double>> m_sol;

    // grille de representation
    unsigned int m_nX; // nombre de cellules selon direction X
    unsigned int m_nY; // idem precedemment
    double m_lX;
    double m_lY;
    std::vector<dealii::Point<2>> m_grid;
    std::vector<std::pair<std::complex<double>,std::complex<double>>> m_data;
};

#endif