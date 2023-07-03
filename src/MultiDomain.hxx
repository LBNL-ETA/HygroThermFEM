#pragma once

#include "Domain.hxx"
#include "ThermalDomain.hxx"
#include "MoistureDomain.hxx"
#include "MaterialMissingProperties.hxx"

namespace Timesteps
{
    class TimestepObserver;
}

namespace HygroThermFEM
{
    struct MultiDomain
    {
        MultiDomain() = default;
        MultiDomain(bool performThermal, bool performMoisture);

        static std::vector<double> property(Variable property);

        //! \brief Creates element with material reference
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        void createElement(size_t index1,
                           size_t index2,
                           size_t index3,
                           size_t index4,
                           const std::string & materialName);

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! @param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

        //! \brief Assign observer to thermal part of the engine.
        //!
        //! @param observer observer that will be notified about each timestep
        void subscribeThermal(Timesteps::TimestepObserver * observer);

        //! \brief Unsubscribe from thermal notifications
        //!
        //! @param observer observer that will be notified about each timestep
        void unsubscribeThermal(Timesteps::TimestepObserver * observer);

        //! \brief Assign observer to moisture part of the engine.
        //!
        //! @param observer observer that will be notified about each timestep
        void subscribeMoisture(Timesteps::TimestepObserver * observer);

        //! \brief Unsubscribe from moisture notifications
        //!
        //! @param observer observer that will be notified about each timestep
        void unsubscribeMoisture(Timesteps::TimestepObserver * observer);

        //! @brief Deletes Geometry and boundary conditions
        void clearModel();

        ThermalDomain thermalDomain;
        MoistureDomain moistureDomain;
        bool simulateThermal{true};
        bool simulateMoisture{true};
    };

}   // namespace HygroThermFEM
