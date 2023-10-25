#include "acoustic_bem.hpp"

AcousticBEM::AcousticBEM(const std::complex<double> _kappa): m_kernel(_kappa),m_kappa(_kappa),m_fe(1),m_dof_handler(m_triangulation),m_mapping(3)
{
    // ecriture maillage ?
    m_write_mesh = true;

    // ordres de quadrature
    m_qorder_singular = 7;
    m_qorder_close    = 6;
    m_qorder_far      = 3;
    m_close_coef      = 2.0;

    // representation
    m_lX = 2.0;
    m_lY = 2.0;
    m_nX = 200;
    m_nY = 200;
    m_grid.clear();
    m_data.clear();
}

void AcousticBEM::run()
{
    setup_grids();
    setup_system();
    assemble_matrices();
    solve();
    radiate();
}

void AcousticBEM::setup_grids()
{
    std::cout << "- Setup grid... " << std::flush;
    dealii::GridGenerator::hyper_sphere<2>(m_triangulation,dealii::Point<2>(0,0),0.5);
    m_triangulation.refine_global(4);
    // m_triangulation.refine_global(8);

    if(m_write_mesh)
    {
        dealii::GridOut m_gout;
        std::ofstream os("../data/mesh.vtk");
        m_gout.write_vtk(m_triangulation,os);
        os.close();
    }
    std::cout << "Done." << std::endl;
}

void AcousticBEM::setup_system()
{
    std::cout << "- Setup system... " << std::flush;
    m_dof_handler.distribute_dofs(m_fe);
    m_lhs.reinit(m_dof_handler.n_dofs(),m_dof_handler.n_dofs());
    m_rhs.reinit(m_dof_handler.n_dofs());
    m_sol.reinit(m_dof_handler.n_dofs());
    std::cout << "Done." << std::endl;
}

void AcousticBEM::assemble_matrices()
{
    std::cout << "- Assemble matrices... " << std::flush;
    dealii::QGauss<1> quadrature_far(m_qorder_far);
    dealii::QGauss<1> quadrature_close(m_qorder_close);
    dealii::QGauss<1> quadrature_singular(m_qorder_singular);
    dealii::FEValues<1,2> tfe_v_far(m_mapping,m_fe,quadrature_far,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    dealii::FEValues<1,2> tfe_v_close(m_mapping,m_fe,quadrature_close,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    dealii::FEValues<1,2> tfe_v_singular(m_mapping,m_fe,quadrature_singular,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    dealii::FEValues<1,2> sfe_v_far(m_mapping,m_fe,quadrature_far,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    dealii::FEValues<1,2> sfe_v_close(m_mapping,m_fe,quadrature_close,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    const unsigned int n_q_far   = tfe_v_far.n_quadrature_points;
    const unsigned int n_q_close = sfe_v_close.n_quadrature_points;

    assert(n_q_far == sfe_v_far.n_quadrature_points);
    assert(n_q_close == sfe_v_close.n_quadrature_points);

    const unsigned int dofs_per_cell = m_fe.n_dofs_per_cell();

    std::vector<dealii::types::global_dof_index> t_local_dof_indices(dofs_per_cell);
    std::vector<dealii::types::global_dof_index> s_local_dof_indices(dofs_per_cell);

    dealii::FullMatrix<std::complex<double>> cell_lhs(dofs_per_cell,dofs_per_cell);
    dealii::Vector<std::complex<double>> cell_rhs(dofs_per_cell);

    // double loop, non-optimized as all the quadrature rules could be precomputed 
    for(const auto &tcell: m_dof_handler.active_cell_iterators())
    {
        tcell->get_dof_indices(t_local_dof_indices);
        for(const auto &scell: m_dof_handler.active_cell_iterators())
        {
            scell->get_dof_indices(s_local_dof_indices);
            //
            auto position = isCloseInteraction(tcell,scell);
            switch(position)
            {
                case 0: // interaction lointaine
                    buildNonSingularCellMatrix(tcell,scell,tfe_v_far,sfe_v_far,cell_lhs);
                    break;
                case 1: // interaction proche
                    buildNonSingularCellMatrix(tcell,scell,tfe_v_close,sfe_v_close,cell_lhs);
                    break;
                case 2: // interaction singuliere
                    buildSingularCellMatrix(tcell,scell,tfe_v_singular,cell_lhs);
                    break;
                default:
                    std::cout << "Unknown cell -> cell interaction: " << position << std::endl;
                    std::exit(EXIT_FAILURE);
            }
            
            for(unsigned int i=0; i<dofs_per_cell; ++i)
            {
                for(unsigned int j=0; j<dofs_per_cell; ++j)
                {
                    m_lhs(t_local_dof_indices[i],s_local_dof_indices[j]) += cell_lhs(i,j);
                }
            }
        }
        // right-hand-side
        buildRHSCellVector(tcell,tfe_v_far,cell_rhs);
        for(unsigned int i=0; i<dofs_per_cell; ++i)
        {
            m_rhs(t_local_dof_indices[i]) += cell_rhs(i);
        }
    }
    std::cout << "Done." << std::endl;
}

void AcousticBEM::buildRHSCellVector(const dealii::Triangulation<1,2>::cell_iterator &_tcell,dealii::FEValues<1,2> &_tfe,dealii::Vector<std::complex<double>> &_cvector)
{
    _cvector = 0;
    //
    _tfe.reinit(_tcell);
    for(unsigned int iq_t: _tfe.quadrature_point_indices())
    {
        auto &xq_t = _tfe.quadrature_point(iq_t);
        auto pw = -1.0*std::exp(std::complex<double>(0.0,1.0)*m_kappa*xq_t(0));
        for(unsigned int i_t:_tfe.dof_indices())
        {
            _cvector(i_t) += _tfe.shape_value(i_t,iq_t)*pw*_tfe.JxW(iq_t);
        }
    }
}

void AcousticBEM::buildSingularCellMatrix(const dealii::Triangulation<1,2>::cell_iterator &_tcell,const dealii::Triangulation<1,2>::cell_iterator &_scell,dealii::FEValues<1,2> &_tfe,dealii::FullMatrix<std::complex<double>> &_cmatrix)
{
    /////////////////////////////////////////
    // integrate the singular interactions //
    /////////////////////////////////////////
    _cmatrix = 0;
    _tfe.reinit(_tcell);
    // loop over all 't'arget quadrature points
    for(unsigned int iq_t: _tfe.quadrature_point_indices())
    {
        //
        const auto &xq_t = _tfe.quadrature_point(iq_t);
        std::cout << std::endl << _tfe.get_quadrature().point(iq_t)(0) << std::endl;
        dealii::QGaussLogR<1> log_quadrature(m_qorder_close,_tfe.get_quadrature().point(iq_t),1.0,true);
        dealii::FEValues<1,2> log_fevalues(m_mapping,m_fe,log_quadrature,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
        // loop over 's'ource points
        for(unsigned int iq_s: log_fevalues.quadrature_point_indices())
        {
            const auto &xq_s = log_fevalues.quadrature_point(iq_s);
            auto R = dealii::Point<2>(xq_t(0)-xq_s(0),xq_t(1)-xq_s(1)).norm();
            if(R < 1e-12)
            {
                std::cout << "Invalid 'R' value" << std::endl;
                std::exit(1);
            }
            auto kernel_value = m_kernel.single_layer(dealii::Point<2>(xq_t(0)-xq_s(0),xq_t(1)-xq_s(1)));
            if(kernel_value.real() != kernel_value.real())
            {
                std::cout << "Invalid kernel value / R = "<< R << std::endl;
                std::cout << "- xq_t = (" << xq_t(0) << "," << xq_t(1) << ")" << std::endl;
                std::cout << "- xq_s = (" << xq_s(0) << "," << xq_s(1) << ")" << std::endl;
                std::exit(1);
            }
            if(kernel_value.imag() != kernel_value.imag())
            {
                std::cout << "Invalid kernel value / R = "<< R << std::endl;
                std::cout << "- xq_t = (" << xq_t(0) << "," << xq_t(1) << ")" << std::endl;
                std::cout << "- xq_s = (" << xq_s(0) << "," << xq_s(1) << ")" << std::endl;
                std::exit(1);
            }
            for(unsigned int i_t: _tfe.dof_indices())
            {
                for(unsigned int i_s: log_fevalues.dof_indices())
                {
                    _cmatrix(i_t,i_s) += kernel_value*_tfe.shape_value(i_t,iq_t)*log_fevalues.shape_value(i_s,iq_s)*_tfe.JxW(iq_t)*log_fevalues.JxW(iq_s); // could be optimized
                }
            }
        }
    }

}

void AcousticBEM::buildNonSingularCellMatrix(const dealii::Triangulation<1,2>::cell_iterator &_tcell,const dealii::Triangulation<1,2>::cell_iterator &_scell,dealii::FEValues<1,2> &_tfe, dealii::FEValues<1,2> &_sfe, dealii::FullMatrix<std::complex<double>> &_cmatrix)
{
    // reinitialize cell matrix
    _cmatrix = 0;
    // initialize the quadrature data on the 't'arget cell and 's'ource cell
    _tfe.reinit(_tcell);
    _sfe.reinit(_scell);
    // loop over target and source quadrature points
    for(unsigned int iq_t: _tfe.quadrature_point_indices())
    {
        auto &xq_t = _tfe.quadrature_point(iq_t);
        for(unsigned int iq_s: _sfe.quadrature_point_indices())
        {
            auto &xq_s = _sfe.quadrature_point(iq_s);
            // compute the value of the kernel
            auto kernel_value = m_kernel.single_layer(dealii::Point<2>(xq_t(0)-xq_s(0),xq_t(1)-xq_s(1)));
            for(unsigned int i_t: _tfe.dof_indices())
            {
                for(unsigned int i_s: _sfe.dof_indices())
                {
                    _cmatrix(i_t,i_s) += kernel_value*_tfe.shape_value(i_t,iq_t)*_sfe.shape_value(i_s,iq_s)*_tfe.JxW(iq_t)*_sfe.JxW(iq_s); // could be optimized
                }
            }
        }
    }
}

void AcousticBEM::solve()
{
    std::cout << "- Solve... " << std::flush;
    try
    {
        auto lhs = dealii2lapack(m_lhs);
        for(unsigned int i=0; i<lhs.size(); ++i)
        {
            std::cout << lhs[i].r << " " << lhs[i].i << std::endl;
        }
        std::exit(1);
        auto rhs = dealii2lapack(m_rhs);
        int N = rhs.size();
        int NRHS = 1;
        int LDA  = N;
        std::vector<int> IPIV(N);
        int LDB  = N;
        int INFO;
        zgesv_(&N,&NRHS,&(lhs[0]),&LDA,&(IPIV[0]),&(rhs[0]),&LDB,&INFO);
        if(INFO != 0)
        {
            std::cout << "\nLapack solver ZGESV failed with error code '" << INFO << "'" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        for(unsigned int i=0; i<rhs.size(); ++i)
        {
            std::cout << rhs[i].r << " " << rhs[i].i << std::endl;
        }
        std::exit(1);
        m_sol = lapack2dealii(rhs);
    }
    catch(const std::exception& e)
    {
        std::cerr << '\n' << e.what() << '\n';
    }
    
    std::cout << "Done." << std::endl;
}

void AcousticBEM::radiate()
{
    std::cout << "- Radiate... " << std::flush;
    // setup grid
    auto dX = m_lX/m_nX;
    auto dY = m_lY/m_nY;
    auto N  = (m_nX+1)*(m_nY+1);
    m_grid = std::vector<dealii::Point<2>>(N);
    m_data = std::vector<std::pair<std::complex<double>,std::complex<double>>>(N,std::make_pair(std::complex<double>(0,0),std::complex<double>(0,0)));
    unsigned int ind=0;
    for(unsigned int i=0; i<m_nX+1;++i)
    {
        for(unsigned int j=0; j<m_nY+1;++j)
        {
            m_grid[ind](0)   = -m_lX/2+i*dX;
            m_grid[ind++](1) = -m_lY/2+j*dY;
        }
    }
    // setup quadratures for the test space
    dealii::QGauss<1> quadrature_far(m_qorder_far);
    dealii::FEValues<1,2> fe_values(m_mapping,m_fe,quadrature_far,dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);
    const unsigned int dofs_per_cell = m_fe.n_dofs_per_cell();
    std::vector<dealii::types::global_dof_index> s_local_dof_indices(dofs_per_cell);
    //
    for(const auto &scell: m_dof_handler.active_cell_iterators()) // loop over elements, not optimized
    {
        //
        fe_values.reinit(scell);
        scell->get_dof_indices(s_local_dof_indices);
        //
        for(unsigned int iq: fe_values.quadrature_point_indices()) // loop over quadrature nodes
        {
            for(auto i: fe_values.dof_indices()) // loop over dof indices
            {
                auto coef = (fe_values.shape_value(i,iq)*m_sol(s_local_dof_indices[i]))*fe_values.JxW(iq);
                for(unsigned int ir=0; ir<m_grid.size(); ++ir)
                {
                    auto &xr = m_grid[ir];
                    auto &xq = fe_values.quadrature_point(iq);
                    auto kernel_value = m_kernel.single_layer(dealii::Point<2>(xr(0)-xq(0),xr(1)-xq(1)));
                        m_data[ir].first += kernel_value*coef;
                }
            }
        }
    }
    //
    for(unsigned int ir=0; ir<N; ++ir)
    {
        m_data[ir].second = m_data[ir].first + std::exp(std::complex<double>(0,1)*m_kappa*m_grid[ir](0));
    }
    // output to vtk file - scattered field
    std::string name("../data/Us_field.vtk");
    std::ofstream us_ofile; us_ofile.open(name);
    if(!us_ofile.is_open())
    {
        std::cout << "\n\nCould not write file: '" << name << "'" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    us_ofile << "# vtk DataFile Version 4.2" << std::endl;
    us_ofile << "vtk output" << std::endl;
    us_ofile << "ASCII" << std::endl;
    us_ofile << "DATASET STRUCTURED_GRID" << std::endl;
    us_ofile << "DIMENSIONS " << m_nX+1 << " " << m_nY+1 << " 1" << std::endl;
    us_ofile << "POINTS " << N << " double" << std::endl;
    for(unsigned int i=0; i<N; ++i)
    {
        us_ofile << m_grid[i](0) << " " << m_grid[i](1) << " 0" << std::endl;
    }
    us_ofile << "POINT_DATA " << N << std::endl;
    us_ofile << "SCALARS real(Us) double 1" << std::endl;
    us_ofile << "LOOKUP_TABLE default" << std::endl;
    for(unsigned int ir=0; ir<N; ++ir)
    {
        us_ofile << m_data[ir].first.real() << std::endl;
    }
    us_ofile << "SCALARS imag(Us) double 1" << std::endl;
    us_ofile << "LOOKUP_TABLE default" << std::endl;
    for(unsigned int ir=0; ir<N; ++ir)
    {
        us_ofile << m_data[ir].first.imag() << std::endl;
    }
    us_ofile.close();
    // output to vtk file - total field
    name = "../data/Ut_field.vtk";
    std::ofstream ut_ofile; ut_ofile.open(name);
    if(!ut_ofile.is_open())
    {
        std::cout << "\n\nCould not write file: '" << name << "'" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    ut_ofile << "# vtk DataFile Version 4.2" << std::endl;
    ut_ofile << "vtk output" << std::endl;
    ut_ofile << "ASCII" << std::endl;
    ut_ofile << "DATASET STRUCTURED_GRID" << std::endl;
    ut_ofile << "DIMENSIONS " << m_nX+1 << " " << m_nY+1 << " 1" << std::endl;
    ut_ofile << "POINTS " << N << " double" << std::endl;
    for(unsigned int i=0; i<N; ++i)
    {
        ut_ofile << m_grid[i](0) << " " << m_grid[i](1) << " 0" << std::endl;
    }
    ut_ofile << "POINT_DATA " << N << std::endl;
    ut_ofile << "SCALARS real(Ut) double 1" << std::endl;
    ut_ofile << "LOOKUP_TABLE default" << std::endl;
    for(unsigned int ir=0; ir<N; ++ir)
    {
        ut_ofile << m_data[ir].second.real() << std::endl;
    }
    ut_ofile << "SCALARS imag(Ut) double 1" << std::endl;
    ut_ofile << "LOOKUP_TABLE default" << std::endl;
    for(unsigned int ir=0; ir<N; ++ir)
    {
        ut_ofile << m_data[ir].second.imag() << std::endl;
    }
    ut_ofile.close();
    std::cout << "Done." << std::endl;
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

std::vector<zlpk> AcousticBEM::dealii2lapack(const dealii::FullMatrix<std::complex<double>> &_mat) const 
{
    auto V = std::vector<zlpk>(_mat.n_cols()*_mat.n_rows());
    for(unsigned int j=0; j<_mat.n_cols(); ++j)
    {
        for(unsigned int i=0; i<_mat.n_rows(); ++i)
        {
            V[j*_mat.n_rows()+i] = _mat(i,j);
        }
    }
    return V;
}

std::vector<zlpk> AcousticBEM::dealii2lapack(const dealii::Vector<std::complex<double>> &_vec) const
{
    auto V = std::vector<zlpk>(_vec.size());
    for(unsigned int i=0; i<V.size(); ++i)
    {
        V[i] = _vec(i);
    }
    return V;
}

dealii::Vector<std::complex<double>> AcousticBEM::lapack2dealii(const std::vector<zlpk> &_vec) const
{
    auto V = dealii::Vector<std::complex<double>>(_vec.size());
    for(unsigned int i=0; i<V.size(); ++i)
    {
        V(i) = _vec[i];
    }
    return V;
}