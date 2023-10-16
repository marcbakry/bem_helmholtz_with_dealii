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
    dealii::GridGenerator::hyper_sphere<2>(m_triangulation);
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

void AcousticBEM::assemble_matrices()
{
    dealii::QGauss<1> quadrature_far(m_qorder_far);
    dealii::QGauss<1> quadrature_close(m_qorder_close);
    dealii::FEValues<1,2> tfe_v_far(m_mapping,m_fe,quadrature_far,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    dealii::FEValues<1,2> tfe_v_close(m_mapping,m_fe,quadrature_close,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    dealii::FEValues<1,2> sfe_v_far(m_mapping,m_fe,quadrature_far,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    dealii::FEValues<1,2> sfe_v_close(m_mapping,m_fe,quadrature_close,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    const unsigned int n_q_far   = tfe_v_far.n_quadrature_points;
    const unsigned int n_q_close = sfe_v_close.n_quadrature_points;

    assert(n_q_far == sfe_v_far.n_quadrature_points);
    assert(n_q_close == sfe_v_close.n_quadrature_points);

    const unsigned int dofs_per_cell = m_fe.n_dofs_per_cell();

    std::vector<dealii::types::global_dof_index> local_of_indices(dofs_per_cell);

    dealii::FullMatrix<std::complex<double>> cell_matrix(dofs_per_cell,dofs_per_cell);
    dealii::Vector<std::complex<double>> cell_rhs(dofs_per_cell);

    // double loop, non-optimized as all the quadrature rules could be precomputed 
    for(const auto &tcell: m_dof_handler.active_cell_iterators())
    {
        for(const auto &scell: m_dof_handler.active_cell_iterators())
        {
            auto position = isCloseInteraction(tcell,scell);
        }
    }
}

void AcousticBEM::solve()
{
    //
}

void AcousticBEM::radiate()
{
    //
}

unsigned int AcousticBEM::isCloseInteraction(const dealii::Triangulation<1,2>::cell_iterator &_cell1,const dealii::Triangulation<1,2>::cell_iterator &_cell2) const
{
    // retourne 'true' si distance(centre_1,centre_2) <= m_close_coef*(longueur_1+longueur_2)
    // longueur <=> longueur de l'arete <=> distance point 1 au point 2
    auto ctr1 = _cell1->center(true);
    auto ctr2 = _cell2->center(true);
    auto diam1 = _cell1->diameter();
    auto diam2 = _cell2->diameter();
    auto R = (ctr2-ctr1).norm();
    if(R > m_close_coef*(diam1+diam2)) // far
    {
        return 0;
    } 
    else if(R > 1e-9) // close non-singular
    {
        return 1;
    }
    else // singular
    {
        return 2;
    }

}