#include <cmath>

#include "SingleDomain.hxx"
#include "FEMunique.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"
#include "NodePool.hxx"
#include "TimestepData.hxx"
#include "MaterialPool.hxx"

namespace HygroThermFEM
{
    SingleDomain::SingleDomain(const BaseVariable property, bool automaticUpdateOfPreviousTimestep) :
        m_Property(property),
        gasCavities(nullptr),
        m_AutomaticUpdatePreviousTimestep(automaticUpdateOfPreviousTimestep)
    {}

    std::vector<NodeFlux> SingleDomain::flux() const
    {
        return m_Elements.flux();
    }

    void SingleDomain::postProcess(std::vector<double> &)
    {
        if(gasCavities == nullptr)
        {
            gasCavities = std::make_unique<EquivalentGasCavities>(m_Elements);
            gasCavities->setGravityVector(m_GravityVector);
        }
        gasCavities->update();
    }

    void SingleDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        m_GravityVector = gravityVector;
        if(gasCavities != nullptr)
        {
            gasCavities->setGravityVector(gravityVector);
        }
    }

}   // namespace HygroThermFEM
