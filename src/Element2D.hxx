#pragma once

#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"

namespace Conrad {

    // Handles linear 2D element (4 nodes)
    class ElementLinear2D {
    public:
        ElementLinear2D(
                Node2D const &t_Node1,
                Node2D const &t_Node2,
                Node2D const &t_Node3,
                Node2D const &t_Node4,
                double const t_Cond,
                double const t_Rho = 1,
                double const t_Cp = 1);

        ElementLinear2D(ElementLinear2D const &t_Element);

        std::vector<size_t> nodeIndexes() const;

        std::vector<std::vector<double> > conductivity() const;

        std::vector<std::vector<double> > RhoCp() const;

    private:
        std::vector<std::vector<double> > calculateConductionMatrix(size_t const t_IntegrationPointIndex) const;

        std::vector<std::vector<double> > calculateRhoCpMatrix(size_t const t_IntegrationPointIndex) const;

        QuadrilateralLinearGlobal2D m_Element;
        QuadrilateralNodes2D m_ElementNodes;
        double m_Cond;
        double m_Rho;
        double m_Cp;
        std::vector<std::vector<double> > m_Conductivity;
        std::vector<std::vector<double> > m_RhoCp;
    };

}