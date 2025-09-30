#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday
{
    // TODO implement the Faraday class, hint - get inspiration from Ampere
    public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}, m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension>& E)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            for(auto xi = m_grid->dual_dom_start(Direction::X)+1; xi < m_grid->dual_dom_end(Direction::X); ++xi )
            {
                B.x(xi) -= m_dt * ( (E.z(xi +1) - E.z(xi - 1)) - (E.y(xi -1) - E.y(xi -1)) ) / (2 * dx) ; 
                B.y(xi) -= m_dt * ( (E.x(xi +1) - E.x(xi - 1)) - (E.z(xi -1) - E.z(xi -1)) ) / (2 * dx) ; 
                B.z(xi) -= m_dt * ( (E.y(xi +1) - E.y(xi - 1)) - (E.x(xi -1) - E.x(xi -1)) ) / (2 * dx) ; 
            }

        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double m_dt;
};

#endif // HYBRIDIR_FARADAY_HPP
