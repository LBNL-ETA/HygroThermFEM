#pragma once

#include <cstddef>

namespace HygroThermFEM
{
    class SimulationProperties
    {
    public:
        static SimulationProperties & Instance();

        //! \brief Sets iteration parameters for use in every Finite Element Method domain defined
        //! in this engine.
        //!
        //! \param relaxationParameter Number that is > 0 and <=1 that will be used in
        //! iterations to recalculate state variable for start in the next iteration.
        //! \param errorTolerance Tolerance for which solution will be acceptable.
        //! \param maxNumberOfIterations Number of iterations for which solver will attempt to reach
        //! solution with given parameters.
        void setIterationParameters(double relaxationParameter,
                                    double errorTolerance,
                                    std::size_t maxNumberOfIterations);

        //! \brief Resets iteration parameters back to default.
        void resetIterationParameters();

        //! \brief Sets simulation parameters for use in thermal and moisture domains. It defines
        //! what parts of equations will be excluded from finite element model calculations.
        //!
        //! \param excludeWaterLiquidTransportation If set to true, it will exclude water liquid
        //! transportation from mass transfer equation.
        //! \param excludeHeatOfEvaporation If set to true, it will exclude heat of evaporation
        //! calculation from heat transfer equation.
        //! \param excludeCapillaryConduction If set to true, it will exclude capillary conduction
        //! calculation from heat transfer equation.
        //! \param excludeVaporDiffusionConduction If set to true, it will exclude water vapor
        //! diffusion conduction from heat transfer equation.
        //! \param thermalConductivityMoistureAndTemperatureDependent Include moisture and
        //! temperature dependency of thermal conductivity
        void setCalculationParameters(bool excludeWaterLiquidTransportation,
                                      bool excludeHeatOfEvaporation,
                                      bool excludeCapillaryConduction,
                                      bool excludeVaporDiffusionConduction,
                                      bool thermalConductivityMoistureAndTemperatureDependent);

        //! \brief Resets calculation parameters back to default.
        void resetCalculationParameters();

        //! \brief Resets all parameters back to default.
        void reset();

        //! \brief Returns current relaxation parameter value.
        double relaxationParamter() const;

        //! \brief Returns current error tolerance.
        double errorToleranceMoisture() const;

        double errorToleranceTemperature() const;

        //! \brief Returns current maximum allowed number of iterations.
        std::size_t maxNumberOfIterations() const;

        //! \brief Defines if liquid transportation should be taken into account during mass
        //! transfer calculations.
        bool excludeWaterLiquidTransportation() const;

        //! \brief Defines if heat of evaporation should be taken into account during heat transfer
        //! calculations.
        bool excludeHeatOfEvaporation() const;

        //! \brief Defines if capillary conduction should be taken into account during heat transfer
        //! calculations.
        bool excludeCapillaryConduction() const;

        //! \brief Defines if vapor diffusion conduction should be taken into account during heat
        //! transfer calculations.
        bool excludeVaporDiffusionConduction() const;

        //!\ brief Defines if moisture and temperature dependent thermal conductivity will be taken
        //! into account.
        bool thermalConductivityTemperatureAndMoistureDependent() const;

    private:
        SimulationProperties() = default;

        //! \brief Keeps default parameters in single place
        struct DefaultProperties
        {
            const double relaxationParameter{1.0};
            const double errorToleranceMoisture{1e-5};
            const double errorToleranceTemperature{1e-3};
            const size_t maxIterations{50};
            const bool excludeWaterLiquidTransportation{false};
            const bool excludeHeatOfEvaporation{false};
            const bool excludeCapillaryConduction{false};
            const bool excludeVaporDiffusionConduction{false};
            const bool thermalConductivityMoistureAndTemperatureDependent{false};
        } defaultProperties;

        double m_RelaxationParameter{defaultProperties.relaxationParameter};
        double m_ErrorToleranceMoisture{defaultProperties.errorToleranceMoisture};
        double m_ErrorToleranceTemperature{defaultProperties.errorToleranceTemperature};
        size_t m_MaxNumberOfIterations{defaultProperties.maxIterations};
        bool m_ExcludeWaterLiquidTransportation{defaultProperties.excludeWaterLiquidTransportation};
        bool m_ExcludeHeatOfEvaporation{defaultProperties.excludeHeatOfEvaporation};
        bool m_ExcludeCapillaryConduction{defaultProperties.excludeCapillaryConduction};
        bool m_ExcludeVaporDiffusionConduction{defaultProperties.excludeVaporDiffusionConduction};
        bool m_ThermalConductivityMoistureAndTemperatureDependent{defaultProperties.thermalConductivityMoistureAndTemperatureDependent};
    };
}   // namespace HygroThermFEM
