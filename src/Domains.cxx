#include <cmath>
#include <algorithm>
#include <iostream>

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
    FenestrationCommon::SquareMatrix IDomain::steadyStateLeftHandSide()
    {
        auto condMat = m_Elements.conductanceMatrix();
        const auto h = m_BCs.HMatrix();
        condMat += h;

        return condMat;
    }

    std::vector<double> IDomain::steadyStateRightHandSide() const
    {
        return m_BCs.RVector();
    }

    FenestrationCommon::SquareMatrix IDomain::transientM_K_H_Matrix(const double t_DTime)
    {
        auto M = m_Elements.getLumpedMass(t_DTime);
        auto M_K_H = m_Elements.conductanceMatrix();
        M_K_H = M_K_H.addDiagonal(M);
        M_K_H += m_BCs.HMatrix();

        return M_K_H;
    }

    std::vector<double>
      IDomain::transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                    const double t_DTime)
    {
        std::vector<double> M{m_Elements.getLumpedMass(t_DTime)};
        auto R = m_BCs.RVector() + m_Elements.RVector();

        auto B = t_PreviousSolution * M + R;

        return B;
    }

    std::vector<double> IDomain::steadyState()
    {
        auto B = steadyStateRightHandSide();
        return CLinearSolver::solveEigen(steadyStateLeftHandSide(), B);
    }

    SingleSolution IDomain::transient(const std::vector<double> & currentStateValues,
                                      const double t_DTime)
    {
        std::vector<double> solution;
        bool converged{false};
        double currentDTime{t_DTime};
        while(!converged)
        {
            auto current = transientTimestep(currentStateValues, currentDTime);
            solution = current.first;
            converged = current.second;
            if(!converged)
            {
                currentDTime = currentDTime / 2.0;
            }
        }

        return {solution, currentDTime};
    }

    std::pair<std::vector<double>, bool>
      IDomain::transientTimestep(const std::vector<double> & currentStateValues,
                                 const double t_DTime)
    {
        auto A = transientM_K_H_Matrix(t_DTime);

        // This is just for debugging purposes since Eigen vector is invisible.
        auto test = A.toVector();

        auto B = transientMT_R_Vector(currentStateValues, t_DTime);

        // CLinearSolver aSolver;

        std::vector<double> solution;
        bool converged{false};
        bool stopIterations{false};

        if(isLinear())
        {
            solution = CLinearSolver::solveEigen(A, B);
            converged = true;
        }
        else
        {
            solution = currentStateValues;

            auto currentNorm = norm(solution);

            size_t numOfIterations = 0;

            while(!stopIterations && !converged)
            {
                double previousNorm = currentNorm;
                auto temp = A * solution;
                temp = B - temp;

                auto dU = CLinearSolver::solveEigen(A, temp);

                solution = solution + dU;

                currentNorm = norm(solution);

                ++numOfIterations;

                m_BCs.updateNodeValues(solution, m_Property, m_AutomaticUpdatePreviousTimestep);

                A = transientM_K_H_Matrix(t_DTime);
                B = transientMT_R_Vector(currentStateValues, t_DTime);

                converged = std::abs(previousNorm - currentNorm) <= ConvergenceError;

                stopIterations = numOfIterations > MaxIterations;
            }
        }

        m_Elements.updateNodeValues(solution, m_Property, m_AutomaticUpdatePreviousTimestep);

        return std::make_pair(solution, converged);
    }

    bool IDomain::isLinear() const
    {
        return m_BCs.isLinear() && m_Elements.isLinear();
    }

    void IDomain::updateNodeValues(const std::vector<double> & values,
        const BaseVariable property, bool updatePreviousValues)
    {
        m_BCs.updateNodeValues(values, property, updatePreviousValues);
        m_Elements.updateNodeValues(values, property, updatePreviousValues);
    }

    IDomain::IDomain(const BaseVariable property, bool automaticUpdateOfPreviousTimestep) :
    m_Property(property), m_AutomaticUpdatePreviousTimestep(automaticUpdateOfPreviousTimestep)
    {}

    std::vector<NodeFlux> IDomain::flux() const
    {
        return m_Elements.flux();
    }

    void IDomain::postProcess(std::vector<double> &) const
    {
        // Default post processing is to do nothing. Inherited classes should add
        // some functionality if necessary.
    }

    void ThermalDomain::createConvectionBC(const size_t index1,
                                           const size_t index2,
                                           double t_ConvectionCoefficient,
                                           double t_AirTemperature)
    {
        m_BCs.assignBC(fem::make_unique<ConvectionBC>(
          index1, index2, t_ConvectionCoefficient, t_AirTemperature));
    }

    void ThermalDomain::createTemperatureBC(const size_t index1,
                                            const size_t index2,
                                            double t_Temp1,
                                            double t_Temp2)
    {
        m_BCs.assignBC(fem::make_unique<TemperatureBC>(index1, index2, t_Temp1, t_Temp2));
    }

    void ThermalDomain::createTemperatureBC(const size_t index1,
                                            const size_t index2,
                                            const double t_Temp)
    {
        m_BCs.assignBC(fem::make_unique<TemperatureBC>(index1, index2, t_Temp));
    }

    void ThermalDomain::createFluxBC(const size_t index1, const size_t index2, const double t_Flux)
    {
        m_BCs.assignBC(fem::make_unique<FluxBC>(index1, index2, t_Flux));
    }

    void ThermalDomain::createBlackBodyRadiationBC(const size_t index1,
                                                   const size_t index2,
                                                   const double t_Emissivity,
                                                   const double t_RadiationTemperature)
    {
        m_BCs.assignBC(fem::make_unique<BlackBodyRadiationBC>(
          index1, index2, t_Emissivity, t_RadiationTemperature));
    }

    void ThermalDomain::createElement(const size_t index1,
                                      const size_t index2,
                                      const size_t index3,
                                      const size_t index4,
                                      const std::string & materialName)
    {
        m_Elements.assignElement(
          fem::make_unique<ElementThermalLinear2D>(index1, index2, index3, index4, materialName));
    }

    ThermalDomain::ThermalDomain(bool automaticUpdatePreviousTimestep) :
    IDomain(BaseVariable::temperature, automaticUpdatePreviousTimestep)
    {}

    void MoistureDomain::createElement(const size_t index1,
                                       const size_t index2,
                                       const size_t index3,
                                       const size_t index4,
                                       const std::string & materialName)
    {
        m_Elements.assignElement(
          fem::make_unique<ElementMoistureLinear2D>(index1, index2, index3, index4, materialName));
    }

    void MoistureDomain::createMoistureBC(const size_t index1,
                                          const size_t index2,
                                          const double t_AirHumidity,
                                          const double t_AirTemperature)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(fem::make_unique<MoisThermFEM::MoistureBC>(
          index1, index2, Material.name(), t_AirHumidity, t_AirTemperature));
    }

    MoistureDomain::MoistureDomain(bool automaticUpdatePreviousTimestep) :
    IDomain(BaseVariable::humidity, automaticUpdatePreviousTimestep)
    {}

    void MoistureDomain::postProcess(std::vector<double> & solution) const
    {
        IDomain::postProcess(solution);
        for(auto & val : solution)
        {
            if(val > 1)
            {
                val = 1;
            }
            if(val < 0)
            {
                val = 0;
            }
        }
    }

}   // namespace MoisThermFEM