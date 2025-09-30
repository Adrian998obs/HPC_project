#ifndef HYBRIDIR_AMPERE_HPP
#define HYBRIDIR_AMPERE_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ampere
{
public:
    Ampere(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension>& J)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            for(auto xi = m_grid->dual_dom_start(Direction::X)+1; xi < m_grid->dual_dom_end(Direction::X); ++xi )
            {
                J.x(xi) = ((B.z(xi-1) - B.z(xi+1)) - (B.y(xi - 1) - B.y(xi + 1))) / (2 * dx); 
                J.y(xi) = ((B.x(xi-1) - B.x(xi+1)) - (B.z(xi - 1) - B.z(xi + 1))) / (2 * dx); 
                J.z(xi) = ((B.y(xi-1) - B.y(xi+1)) - (B.x(xi - 1) - B.x(xi + 1))) / (2 * dx); 
            }

        }
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP
