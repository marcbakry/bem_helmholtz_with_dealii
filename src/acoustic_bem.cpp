#include "acoustic_bem.hpp"

AcousticBEM::AcousticBEM(const std::complex<double> _kappa): m_kernel(_kappa),m_kappa(_kappa),m_fe(1),m_dof_handler(m_triangulation),m_mapping(3)
{
    // ecriture maillage ?
    m_write_mesh = true;

    // ordres de quadrature
    m_qorder_singular = 6;
    m_qorder_close    = 6;
    m_qorder_far      = 3;
    m_close_coef      = 2.0;
}

void AcousticBEM::run()
{
    setup_grids();
}

void AcousticBEM::setup_grids()
{
    dealii::GridGenerator::hyper_sphere(m_triangulation);
    m_triangulation.refine_global(8);

    if(m_write_mesh)
    {
        dealii::GridOut m_gout;
        std::ofstream os("../data/mesh.vtk");
        m_gout.write_vtk(m_triangulation,os);
        os.close();
    }
}

void AcousticBEM::setup_system()
{
    m_dof_handler.distribute_dofs(m_fe);
    m_lhs.reinit(m_dof_handler.n_dofs(),m_dof_handler.n_dofs());
    m_rhs.reinit(m_dof_handler.n_dofs());
    m_sol.reinit(m_dof_handler.n_dofs());
}

bool AcousticBEM::isCloseInteraction(const dealii::Triangulation<1,2>::cell_iterator &_cell1,const dealii::Triangulation<1,2>::cell_iterator &_cell2) const
{
    // retourne 'true' si distance(centre_1,centre_2) <= m_close_coef*(longueur_1+longueur_2)
    // longueur <=> longueur de l'arete <=> distance point 1 au point 2
}