#include <iostream>
#include <vector>
#include "highfive/highfive.hpp"
#include "field.hpp"
#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "ampere.hpp"
#include "faraday.hpp"
#include "boundary_condition.hpp"
#include <cmath>

void test_ampere()
{   
    std::cout << "Running ampère test...\n";
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 1;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);

    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        auto x_y = layout->coordinate(Direction::X, Quantity::By, ix);
        auto x_z = layout->coordinate(Direction::X, Quantity::Bz, ix);
        B.y(ix) = cos(x_y);
        B.z(ix) = sin(x_z);
    }
    boundary_condition->fill(B);

    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    
    Ampere<dimension> ampere{layout};
    ampere(B, J);
    boundary_condition->fill(J);
    {
        std::string filename = "ampere.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Jx", J.x.data());
        file.createDataSet("/Jy", J.y.data());
        file.createDataSet("/Jz", J.z.data());
    }

    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        auto x_y = layout->coordinate(Direction::X, Quantity::By, ix);
        auto x_z = layout->coordinate(Direction::X, Quantity::Bz, ix);
        double delta_Jy = std::abs(J.y(ix) + cos(x_y));
        double delta_Jz = std::abs(J.z(ix) + sin(x_z));
        if (delta_Jy > 1e-2 || delta_Jz > 1e-2)
        {
            std::cout << "Ampere test failed at ix=" << ix << " x=" << x_y << " Jy=" << J.y(ix)
                      << " Jz=" << J.z(ix) << " delta_Jy=" << delta_Jy
                      << " delta_Jz=" << delta_Jz << "\n";
            return;
        }
    }
}

void test_faraday()
{
    std::cout << "Running faraday test...\n";
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 0;
    auto constexpr nppc                          = 100;
    double dt                                    = 0.001;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);


    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Bx, ix);
        B.x(ix) = 3 * x;   
    }
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, ix);
        B.y(ix) = 0.0;
        B.z(ix) = 1.0; 
        E.y(ix) = sin(x);
        E.z(ix) = cos(x);
    }

    Faraday<dimension> faraday{layout, dt};
    faraday(B, Bnew, E);
    boundary_condition->fill(Bnew);
    boundary_condition->fill(E);

    {
        std::string filename = "faraday.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Bnewx", Bnew.x.data());
        file.createDataSet("/Bnewy", Bnew.y.data());
        file.createDataSet("/Bnewz", Bnew.z.data());
    }

    for (auto ix = layout->primal_dom_start(Direction::X); ix <= layout->primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::Bx, ix);
        double delta_Bx = std::abs(Bnew.x(ix) - 3 * x);
        if (delta_Bx > 1e-2)
        {
            std::cout << "Faraday test failed at ix=" << ix << " x=" << x << " Bx=" << Bnew.x(ix)
                      << " delta_Bx=" << delta_Bx << "\n";
            return;
        }
    }
    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, ix);
        double expected_By = - dt * sin(x);
        double expected_Bz = 1.0 - dt * cos(x);
        double delta_By    = std::abs(Bnew.y(ix) - expected_By);
        double delta_Bz    = std::abs(Bnew.z(ix) - expected_Bz);
        if (delta_By > 1e-2 || delta_Bz > 1e-2)
        {
            std::cout << "Faraday test failed at ix=" << ix << " x=" << x << " By=" << Bnew.y(ix)
                      << " Bz=" << Bnew.z(ix) << " delta_By=" << delta_By
                      << " delta_Bz=" << delta_Bz << "\n";
            return;
        }
    }
}

int main()
{
    test_ampere();
    test_faraday();
    return 0;
}