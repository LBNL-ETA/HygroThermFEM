#include <cmath>
#include <algorithm>

#include "Domains.hxx"
#include "FEMunique.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"

using FenestrationCommon::CLinearSolver;

namespace MoisThermFEM
{

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

    std::vector<double> Domain::transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
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

    std::vector<double> Domain::transient(const std::vector<double> & currentStateValues,
                                          const double t_DTime)
    {
        auto A = transientM_K_H_Matrix(t_DTime);
        auto a = A.toVector();
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

                //std::transform(
                //  dU.begin(), dU.end(), solution.begin(), solution.begin(), std::plus<double>());

                solution = solution + dU;

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

	IElementLinear2D * Domain::findElement(const Node2D & t_Node1, const Node2D & t_Node2)
	{
		return m_Elements.findElement(t_Node1, t_Node2);
	}

	Domain::Domain( const Property property ) : m_Property(property) {

	}

	void ThermalDomain::createConvectionBC(const Node2D & t_Node1,
                                    const Node2D & t_Node2,
                                    const double t_ConvectionCoefficient,
                                    const double t_AirTemperature)
    {
    	m_BCs.assignBC(fem::make_unique<ConvectionBC>(
			t_Node1, t_Node2, t_ConvectionCoefficient, t_AirTemperature));
    }

    void ThermalDomain::createTemperatureBC(Node2D & t_Node1,
                                     Node2D & t_Node2,
                                     const double t_Temp1,
                                     const double t_Temp2)
    {
    	m_BCs.assignBC(fem::make_unique<TemperatureBC>(t_Node1, t_Node2, t_Temp1, t_Temp2));
    }

    void ThermalDomain::createTemperatureBC(Node2D & t_Node1, Node2D & t_Node2, const double t_Temp)
    {
    	m_BCs.assignBC(fem::make_unique<TemperatureBC>(t_Node1, t_Node2, t_Temp));
    }

    void ThermalDomain::createFluxBC(Node2D & t_Node1, Node2D & t_Node2, const double t_Flux)
    {
    	m_BCs.assignBC(fem::make_unique<FluxBC>(t_Node1, t_Node2, t_Flux));
    }

    void ThermalDomain::createBlackBodyRadiationBC(const Node2D & t_Node1,
                                            const Node2D & t_Node2,
                                            const double t_Emissivity,
                                            const double t_RadiationTemperature)
    {
    	m_BCs.assignBC(fem::make_unique<BlackBodyRadiationBC>(
			t_Node1, t_Node2, t_Emissivity, t_RadiationTemperature));
    }

    void ThermalDomain::createElement(const Node2D & t_Node1,
                                      const Node2D & t_Node2,
                                      const Node2D & t_Node3,
                                      const Node2D & t_Node4,
                                      const Material & mat)
    {
        m_Elements.assignElement(
          fem::make_unique<ElementThermalLinear2D>(t_Node1, t_Node2, t_Node3, t_Node4, mat));
    }

	ThermalDomain::ThermalDomain() : Domain(Property::temperature) {

	}

	void MoistureDomain::createElement(const Node2D & t_Node1,
                                       const Node2D & t_Node2,
                                       const Node2D & t_Node3,
                                       const Node2D & t_Node4,
                                       const Material & mat)
    {
        m_Elements.assignElement(
          fem::make_unique<ElementMoistureLinear2D>(t_Node1, t_Node2, t_Node3, t_Node4, mat));
    }

	void MoistureDomain::createMoistureBC(const Node2D & t_Node1,
										  const Node2D & t_Node2,
										  const double t_ConvectiveCoefficient,
										  const double t_AirHumidity,
										  const double t_AirTemperature)
	{
		/// Need to pull material for current moisture boundary condition
		auto & Material = m_Elements.findElement(t_Node1, t_Node2)->getMaterial();
		m_BCs.assignBC(fem::make_unique<MoisThermFEM::MoistureBC>(
			t_Node1, t_Node2, t_ConvectiveCoefficient, Material, t_AirHumidity, t_AirTemperature));
	}

	MoistureDomain::MoistureDomain() : Domain(Property::humidity) {

	}

}   // namespace MoisThermFEM