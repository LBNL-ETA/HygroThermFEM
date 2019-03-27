#include <cmath>

#include "Domains.hxx"
#include "FEMunique.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"

namespace HygroThermFEM
{
    SquareMatrix IDomain::steadyStateLeftHandSide()
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

    SquareMatrix IDomain::transientM_K_H_Matrix(const double t_DTime)
    {
        const auto M = m_Elements.getLumpedMass(t_DTime);
        auto M_K_H = m_Elements.conductanceMatrix();
        M_K_H = M_K_H.addDiagonal(M);
        M_K_H += m_BCs.HMatrix();

        return M_K_H;
    }

    std::vector<double>
      IDomain::transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                    const double t_DTime)
    {
        const std::vector<double> M{m_Elements.getLumpedMass(t_DTime)};
        const auto R = m_BCs.RVector() + m_Elements.RVector();

        auto B = t_PreviousSolution * M + R;

        return B;
    }

    std::vector<double> IDomain::steadyState()
    {
        const auto B = steadyStateRightHandSide();
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
            const auto current = transientTimestep(currentStateValues, currentDTime);
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
        const auto RelaxParameter = 1.0;
        auto A = transientM_K_H_Matrix(t_DTime);

        // This is just for debugging purposes.
        // auto testA = A.toVector();

        auto B = transientMT_R_Vector(currentStateValues, t_DTime);

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
                const double previousNorm = currentNorm;
                auto temp = A * solution;
                temp = B - temp;

                /*
                auto testA = A.toVector();

                std::cout.precision(18);
                std::cout << std::endl << "------------------------------------------";
                std::cout << std::endl;
                for(const auto & row : testA) {
                    for(const auto & val : row) {
                        std::cout << val << ",";
                    }
                    std::cout << std::endl;
                }

                std::cout << std::endl << std::endl;
                std::cout << std::endl << "------------------------------------------" << std::endl;

                for(const auto & val: B) {
                    std::cout << val << std::endl;
                }

                std::cout << std::endl;
                std::cout << std::endl << "------------------------------------------";
                */

                auto dU = CLinearSolver::solveEigen(A, temp);

                solution = solution + dU * RelaxParameter;

                currentNorm = norm(solution);

                ++numOfIterations;

                m_BCs.updateNodeValues(solution, m_Property, m_AutomaticUpdatePreviousTimestep);
                m_Elements.updateNodeValues(
                  solution, m_Property, m_AutomaticUpdatePreviousTimestep);

                A = transientM_K_H_Matrix(t_DTime);
                // test = A.toVector();
                B = transientMT_R_Vector(currentStateValues, t_DTime);

                converged = (std::abs(previousNorm - currentNorm) / (currentNorm + 1e-6))
                            <= (ConvergenceError * RelaxParameter);

                stopIterations = numOfIterations > (MaxIterations / RelaxParameter);

                postProcess(solution);
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
                                   const BaseVariable property,
                                   bool updatePreviousValues)
    {
        m_BCs.updateNodeValues(values, property, updatePreviousValues);
        m_Elements.updateNodeValues(values, property, updatePreviousValues);
    }

    IDomain::IDomain(const BaseVariable property, bool automaticUpdateOfPreviousTimestep) :
        m_Property(property),
        m_AutomaticUpdatePreviousTimestep(automaticUpdateOfPreviousTimestep)
    {}

    std::vector<NodeFlux> IDomain::flux() const
    {
        return m_Elements.flux();
    }

    void IDomain::postProcess(std::vector<double> &)
    {
        // Default post processing is to do nothing. Inherited classes should add
        // some functionality if necessary.
    }

    void ThermalDomain::createConvectionBCFixedHc(const size_t index1,
                                                  const size_t index2,
                                                  const double t_AirTemperature,
                                                  const double t_ConvectionCoefficient,
                                                  const double t_AirHumidity,
                                                  const bool t_CalculateMoisture)
    {
        m_BCs.assignBC(fem::make_unique<ConstantConvectionBC>(index1,
                                                              index2,
                                                              t_AirTemperature,
                                                              t_ConvectionCoefficient,
                                                              t_AirHumidity,
                                                              t_CalculateMoisture));
    }

    void ThermalDomain::createConvectionBCVariableHc(const size_t index1,
                                                     const size_t index2,
                                                     const double t_AirTemperature,
                                                     const double t_AirHumidity,
                                                     const bool t_CalculateMoisture)
    {
        m_BCs.assignBC(fem::make_unique<VariableConvectionBC>(
          index1, index2, t_AirTemperature, t_AirHumidity, t_CalculateMoisture));
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

    void ThermalDomain::createSimplifiedRadiationBC(const size_t index1,
                                                    const size_t index2,
                                                    const double t_RadiationCoefficient,
                                                    const double t_RadiationTemperature)
    {
        m_BCs.assignBC(fem::make_unique<SimplifiedRadiationBC>(
          index1, index2, t_RadiationCoefficient, t_RadiationTemperature));
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

    void ThermalDomain::postProcess(std::vector<double> & solution)
    {
        IDomain::postProcess(solution);
        if(frameCavities == nullptr)
        {
            frameCavities =
              std::unique_ptr<EquivalentFrameCavities>(new EquivalentFrameCavities(m_Elements));
        }
        frameCavities->update();
    }

    ThermalDomain::ThermalDomain(bool automaticUpdatePreviousTimestep) :
        IDomain(BaseVariable::temperature, automaticUpdatePreviousTimestep),
        frameCavities(nullptr)
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

    void MoistureDomain::createMoistureBCVariableHc(const size_t index1,
                                                    const size_t index2,
                                                    const double t_AirHumidity,
                                                    const double t_AirTemperature)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(fem::make_unique<HygroThermFEM::MoistureBCVariableHc>(
          index1, index2, Material.name(), t_AirHumidity, t_AirTemperature));
    }

    void MoistureDomain::createMoistureBCFixedHc(const size_t index1,
                                                 const size_t index2,
                                                 const double t_AirTemperature,
                                                 const double t_ConvectiveFilmCoefficient,
                                                 const double t_AirHumidity)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(
          fem::make_unique<HygroThermFEM::MoistureBCFixedHc>(index1,
                                                             index2,
                                                             Material.name(),
                                                             t_AirHumidity,
                                                             t_AirTemperature,
                                                             t_ConvectiveFilmCoefficient));
    }

    MoistureDomain::MoistureDomain(bool automaticUpdatePreviousTimestep) :
        IDomain(BaseVariable::humidity, automaticUpdatePreviousTimestep)
    {}

    void MoistureDomain::postProcess(std::vector<double> & solution)
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

}   // namespace HygroThermFEM
