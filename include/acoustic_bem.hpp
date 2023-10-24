#ifndef ACOUSTIC_BEM_HPP
#define ACOUSTIC_BEM_HPP

#include <fstream>
#include <cassert>
#include <complex>
#include <cmath>

#include "defs_dealii.hpp"
#include "helmholtz_kernel.hpp"

#include "cblas.h" // use LAPACK

extern "C"
{
    typedef struct{double r,i;} __CLPK_doublecomplex;
}
class zlpk : public __CLPK_doublecomplex
{
public:
    zlpk() {r=0; i=0;};
    zlpk(double v, double w=0) {r=v; i=w;};
    zlpk(std::complex<double> v) {r=v.real(); i=v.imag();};
    operator std::complex<double>() const {return std::complex<double>(r,i);};
};

extern "C"
{
    void zgesv_(int* N, int* NRHS, zlpk* A, int* LDA, int* IPIV, zlpk* B, int* LDB, int* INFO);
}

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

    std::vector<zlpk> dealii2lapack(const dealii::FullMatrix<std::complex<double>> &_mat) const;
    std::vector<zlpk> dealii2lapack(const dealii::Vector<std::complex<double>> &_vec) const;
    dealii::Vector<std::complex<double>> lapack2dealii(const std::vector<zlpk> &_vec) const;

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

    // grille de representation
    unsigned int m_nX; // nombre de cellules selon direction X
    unsigned int m_nY; // idem precedemment
    double m_lX;
    double m_lY;
    std::vector<dealii::Point<2>> m_grid;
    std::vector<std::complex<double>> m_data;
};

#endif