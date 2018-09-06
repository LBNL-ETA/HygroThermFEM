#include <cmath>
#include <algorithm>

#include "Domains.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "VectorOperators.hxx"

using FenestrationCommon::CLinearSolver;

namespace MoisThermFEM
{
    Domain::Domain(const Property property) : m_Property{property}
    {}

    FenestrationCommon::SquareMatrix Domain::steadyStateLeftHandSide()
    {
        auto condMat = m_Elements.conductanceMatrix();
        const auto h = m_BCs.HMatrix();
        condMat = condMat + h;

        return condMat;
    }

    std::vector<double> Domain::steadyStateRightHandSide() const
    {
        return m_BCs.RVector();
    }

    FenestrationCommon::SquareMatrix Domain::transientM_K_H_Matrix(const double t_DTime)
    {
        auto M = m_Elements.getLumpedMass(t_DTime);
        // auto M = m_Elements.getMassMatrix( t_DTime );
        auto M_K_H = m_Elements.conductanceMatrix();
        M_K_H = M_K_H.addDiagonal(M);
        // M_K_H += M;
        M_K_H += m_BCs.HMatrix();

        return M_K_H;
    }

    std::vector<double> Domain::transientMT_R_Vector(std::vector<double> & t_PreviousSolution,
                                                     const double t_DTime)
    {
        std::vector<double> M{m_Elements.getLumpedMass(t_DTime)};
        auto R = m_BCs.RVector();

        auto B = t_PreviousSolution * M + R;

        return B;
    }

    std::vector<double> Domain::steadyState()
    {
        auto B = steadyStateRightHandSide();
        return CLinearSolver::solveEigen(steadyStateLeftHandSide(), B);
    }

    std::vector<double> Domain::transient(std::vector<double> & currentStateValues,
                                          const double t_DTime)
    {
        auto A = transientM_K_H_Matrix(t_DTime);
        auto B = transientMT_R_Vector(currentStateValues, t_DTime);

        // CLinearSolver aSolver;

        std::vector<double> solution;

        if(isLinear())
        {
            solution = CLinearSolver::solveEigen(A, B);
        }
        else
        {
            solution = currentStateValues;

            auto error = std::numeric_limits<double>::max();

            size_t numOfIterations = 0;

            while(error > ConvergenceError)
            {
                auto temp = A * solution;
                temp = B - temp;

                /// Seems that DH can be avoided. Same solution is achieved faster without it. Topaz
                /// does have this implementation. Will keep it commented in case we want to test it
                /// in future when new kind of boundary conditions are introduced (Simon) auto DH =
                /// transientDH_Matrix( ); DH = A + DH;

                /// auto dU = aSolver.solveSystem( DH, temp );

                auto dU = CLinearSolver::solveEigen(A, temp);

                error = norm(dU);

                // for( auto i = 0u; i < solution.size(); ++i ) {
                // 	solution[i] += dU[i];
                // }

                std::transform(
                  dU.begin(), dU.end(), solution.begin(), solution.begin(), std::plus<double>());

                ++numOfIterations;

                m_BCs.updateNodeValues(solution, m_Property);

                A = transientM_K_H_Matrix(t_DTime);
                B = transientMT_R_Vector(currentStateValues, t_DTime);

                if(numOfIterations > MaxIterations)
                {
                    throw std::runtime_error("Solution failed to converge.");
                }
            }
        }

        return solution;
    }

    bool Domain::isLinear() const
    {
        return m_BCs.isLinear() && m_Elements.isLinear();
    }

    void Domain::updateNodeValues(const std::vector<double> & values, const Property property)
    {
        m_BCs.updateNodeValues(values, property);
        m_Elements.updateNodeValues(values, property);
    }

    void Domain::createConvectionBC(const Node2D & t_Node1,
                                    const Node2D & t_Node2,
                                    const double t_ConvectionCoefficient,
                                    const double t_AirTemperature)
    {
        m_BCs.createConvectionBC(t_Node1, t_Node2, t_ConvectionCoefficient, t_AirTemperature);
    }

    void Domain::createTemperatureBC(Node2D & t_Node1,
                                     Node2D & t_Node2,
                                     const double t_Temp1,
                                     const double t_Temp2)
    {
        m_BCs.createTemperatureBC(t_Node1, t_Node2, t_Temp1, t_Temp2);
    }

    void Domain::createTemperatureBC(Node2D & t_Node1, Node2D & t_Node2, const double t_Temp)
    {
        m_BCs.createTemperatureBC(t_Node1, t_Node2, t_Temp);
    }

    void Domain::createFluxBC(Node2D & t_Node1, Node2D & t_Node2, const double t_Flux)
    {
        m_BCs.createFluxBC(t_Node1, t_Node2, t_Flux);
    }

    void Domain::createBlackBodyRadiationBC(const Node2D & t_Node1,
                                            const Node2D & t_Node2,
                                            const double t_Emissivity,
                                            const double t_RadiationTemperature)
    {
        m_BCs.createBlackBodyRadiationBC(t_Node1, t_Node2, t_Emissivity, t_RadiationTemperature);
    }

    void Domain::createMoistureBC(const Node2D & t_Node1,
                                  const Node2D & t_Node2,
                                  const double t_ConvectiveCoefficient,
                                  const double t_AirHumidity,
                                  const double t_AirTemperature)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = findElement(t_Node1, t_Node2)->getMaterial();
        m_BCs.createMoistureBC(
          t_Node1, t_Node2, t_ConvectiveCoefficient, Material, t_AirHumidity, t_AirTemperature);
    }

    void Domain::createThermalElement(const Node2D & t_Node1,
                                      const Node2D & t_Node2,
                                      const Node2D & t_Node3,
                                      const Node2D & t_Node4,
                                      const Material & mat)
    {
        m_Elements.createThermalElement(t_Node1, t_Node2, t_Node3, t_Node4, mat);
    }

    void Domain::createMoistureElement(const Node2D & t_Node1,
                                       const Node2D & t_Node2,
                                       const Node2D & t_Node3,
                                       const Node2D & t_Node4,
                                       const Material & mat)
    {
        m_Elements.createMoistureElement(t_Node1, t_Node2, t_Node3, t_Node4, mat);
    }

    IElementLinear2D * Domain::findElement(const Node2D & t_Node1, const Node2D & t_Node2)
    {
        return m_Elements.findElement(t_Node1, t_Node2);
    }

    /// FenestrationCommon::SquareMatrix< double > Domain::transientDH_Matrix() {
    /// 	return m_BCs.DHMatrix();;
    /// }
}   // namespace MoisThermFEM