#include <cmath>

#include "Domains.hxx"
#include "FEMunique.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"
#include "NodePool.hxx"

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

    SquareMatrix IDomain::transientM_K_H_Matrix(const double t_DTime, const size_t timestepIndex)
    {
        const auto M = m_Elements.getLumpedMass(t_DTime);
        auto M_K_H = m_Elements.conductanceMatrix();
        M_K_H = M_K_H.addDiagonal(M);
        M_K_H += m_BCs.HMatrix(timestepIndex);

        return M_K_H;
    }

    std::vector<double>
      IDomain::transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                    const double t_DTime,
                                    const size_t timestepIndex)
    {
        const std::vector<double> M{m_Elements.getLumpedMass(t_DTime)};
        const auto R = m_BCs.RVector(timestepIndex) + m_Elements.RVector();

        auto B = t_PreviousSolution * M + R;

        return B;
    }

    std::vector<double> IDomain::steadyState()
    {
        const auto B = steadyStateRightHandSide();
        const auto A = steadyStateLeftHandSide();
        // const auto test = A.toVector();
        return CLinearSolver::solveEigen(A, B);
    }

    SingleSolution IDomain::transient(const std::vector<double> & currentStateValues,
                                      const double t_DTime,
                                      const size_t timestepIndex)
    {
        std::vector<double> solution;
        bool converged{false};
        size_t maxDivisions{3u};
        size_t currentDivision{0u};
        double currentDTime{t_DTime};
        double totalTime{0};
        auto stateVariables{currentStateValues};
        // In case program failed to converge, it will cut down step to smaller one and will perform
        // multiple consecutive simulations in order to achieve solution at requested timestep.
        while(totalTime < t_DTime)
        {
            const auto current = transientTimestep(stateVariables, currentDTime, timestepIndex);
            solution = current.first;
            converged = current.second;
            if(!converged)
            {
                currentDTime = currentDTime / 10;
                ++currentDivision;
                if(currentDivision > maxDivisions)
                {
                    throw std::runtime_error("Solution failed to converge.");
                }
            }
            else
            {
                stateVariables = solution;
                totalTime += currentDTime;
            }
        }

        return {solution, t_DTime};
    }

    std::pair<std::vector<double>, bool>
      IDomain::transientTimestep(const std::vector<double> & currentStateValues,
                                 const double t_DTime,
                                 const size_t timestepIndex)
    {
        const auto RelaxParameter = SimulationProperties::Instance().relaxationParamter();
        const auto ConvergenceError = SimulationProperties::Instance().errorTolerance();
        const auto MaxIterations = SimulationProperties::Instance().maxNumberOfIterations();

        auto A = transientM_K_H_Matrix(t_DTime, timestepIndex);

        auto B = transientMT_R_Vector(currentStateValues, t_DTime, timestepIndex);

        std::vector<double> solution;
        bool converged{false};
        bool stopIterations{false};

        if(isLinear())
        {
            solution = CLinearSolver::solveEigen(A, B);
            postProcess(solution);
            converged = true;
        }
        else
        {
            solution = currentStateValues;
            std::vector<double> normSolution{currentStateValues};

            auto currentNorm = norm(solution);

            size_t numOfIterations = 0;

            while(!stopIterations && !converged)
            {
                const double previousNorm = currentNorm;
                auto temp = A * solution;
                temp = B - temp;

                auto dU = CLinearSolver::solveEigen(A, temp);

                solution = solution + dU * RelaxParameter;
                normSolution = solution + dU;

                postProcess(solution);
                postProcess(normSolution);

                currentNorm = norm(normSolution);

                ++numOfIterations;

                NodePool::Instance().updateNodeValues(
                  solution, m_Property, m_AutomaticUpdatePreviousTimestep);

                A = transientM_K_H_Matrix(t_DTime, timestepIndex);
                B = transientMT_R_Vector(currentStateValues, t_DTime, timestepIndex);

                converged = (std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12))
                            <= ConvergenceError;

                stopIterations = numOfIterations > MaxIterations;
            }
        }

        NodePool::Instance().updateNodeValues(
          solution, m_Property, m_AutomaticUpdatePreviousTimestep);

        return std::make_pair(solution, converged);
    }

    bool IDomain::isLinear() const
    {
        return m_BCs.isLinear() && m_Elements.isLinear();
    }

    IDomain::IDomain(const BaseVariable property, bool automaticUpdateOfPreviousTimestep) :
        m_Property(property), m_AutomaticUpdatePreviousTimestep(automaticUpdateOfPreviousTimestep)
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

    void
      ThermalDomain::createConvectionBCFixedHc(const size_t index1,
                                               const size_t index2,
                                               const FixedBCHCCoefficients & fixedBCHCCoefficients,
                                               const bool t_CalculateMoisture)
    {
        m_BCs.assignBC(fem::make_unique<ConstantConvectionBC>(
          index1, index2, fixedBCHCCoefficients, t_CalculateMoisture));
    }

    void ThermalDomain::createConvectionBCFixedHc(
      size_t index1,
      size_t index2,
      const std::vector<FixedBCHCCoefficients> & fixedBCHCCoefficients,
      bool calculateMoisture)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(
          fixedBCHCCoefficients.begin(), fixedBCHCCoefficients.end(), [&](const auto & bc) {
              timestepBCs.push_back(
                std::make_unique<ConstantConvectionBC>(index1, index2, bc, calculateMoisture));
          });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createConvectionBCVariableHc(const size_t index1,
                                                     const size_t index2,
                                                     const VariableBCHCCoefficients & varHCCoeff,
                                                     const bool t_CalculateMoisture)
    {
        m_BCs.assignBC(
          std::make_unique<VariableConvectionBC>(index1, index2, varHCCoeff, t_CalculateMoisture));
    }

    void ThermalDomain::createConvectionBCVariableHc(
      size_t index1,
      size_t index2,
      const std::vector<VariableBCHCCoefficients> & varHCCoeff,
      const bool t_CalculateMoisture)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(varHCCoeff.begin(), varHCCoeff.end(), [&](const auto & bc) {
            timestepBCs.push_back(
              std::make_unique<VariableConvectionBC>(index1, index2, bc, t_CalculateMoisture));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createTemperatureBC(const size_t index1,
                                            const size_t index2,
                                            double t_Temp1,
                                            double t_Temp2)
    {
        m_BCs.assignBC(std::make_unique<TemperatureBC>(index1, index2, t_Temp1, t_Temp2));
    }

    void ThermalDomain::createTemperatureBC(size_t index1,
                                            size_t index2,
                                            const std::vector<ConstantBCTemperatures> & temp)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(temp.begin(), temp.end(), [&](const auto & bc) {
            timestepBCs.push_back(
              std::make_unique<TemperatureBC>(index1, index2, bc.Temperature1, bc.Temperature2));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createTemperatureBC(const size_t index1,
                                            const size_t index2,
                                            const double t_Temp)
    {
        m_BCs.assignBC(fem::make_unique<TemperatureBC>(index1, index2, t_Temp));
    }

    void ThermalDomain::createTemperatureBC(size_t index1, size_t index2, std::vector<double> temp)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(temp.begin(), temp.end(), [&](const auto & temperature) {
            timestepBCs.push_back(std::make_unique<TemperatureBC>(index1, index2, temperature));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
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

    void ThermalDomain::createBlackBodyRadiationBC(
      size_t index1, size_t index2, const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(radCoeffs.begin(), radCoeffs.end(), [&](const auto & bc) {
            timestepBCs.push_back(std::make_unique<BlackBodyRadiationBC>(
              index1, index2, bc.Emissivity, bc.Temperature));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createSimplifiedRadiationBC(
      const size_t index1,
      const size_t index2,
      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        m_BCs.assignBC(fem::make_unique<SimplifiedRadiationBC>(index1, index2, linearRadBC));
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
        for(auto & val : solution)
        {
            if(val < Constants::ABSOLUTEZERO)
            {
                val = Constants::ABSOLUTEZERO + 1e-6;
            }
            if(val > 1000)
            {
                val = 1000;
            }
        }
        if(frameCavities == nullptr)
        {
            frameCavities = std::make_unique<EquivalentFrameCavities>(m_Elements);
        }
        frameCavities->update();
    }

    ThermalDomain::ThermalDomain(bool automaticUpdatePreviousTimestep) :
        IDomain(BaseVariable::temperature, automaticUpdatePreviousTimestep), frameCavities(nullptr)
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
                                                    const VariableBCHCCoefficients & varHCCoeff)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(
          fem::make_unique<MoistureBCVariableHc>(index1, index2, Material.name(), varHCCoeff));
    }

    void MoistureDomain::createMoistureBCVariableHc(
      size_t index1, size_t index2, const std::vector<VariableBCHCCoefficients> & varCoeff)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        std::for_each(varCoeff.begin(), varCoeff.end(), [&](const auto & bc) {
            timestepBCs.push_back(
              std::make_unique<MoistureBCVariableHc>(index1, index2, Material.name(), bc));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void MoistureDomain::createMoistureBCFixedHc(
      const size_t index1, const size_t index2, const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(std::make_unique<MoistureBCFixedHc>(
          index1, index2, Material.name(), fixedBchcCoefficients));
    }

    void MoistureDomain::createMoistureBCFixedHc(
      size_t index1,
      size_t index2,
      const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        std::for_each(
          fixedBchcCoefficients.begin(), fixedBchcCoefficients.end(), [&](const auto & bc) {
              timestepBCs.push_back(
                std::make_unique<MoistureBCFixedHc>(index1, index2, Material.name(), bc));
          });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
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
