#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // FixedBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    FixedBCHCCoefficients::FixedBCHCCoefficients(const double airTemperature,
                                                 const double convectionCoefficient) :
        AirTemperature(airTemperature), ConvectionCoefficient(convectionCoefficient)
    {}

    FixedBCHCCoefficients::FixedBCHCCoefficients(const double airTemperature,
                                                 const double convectionCoefficient,
                                                 const double airHumidity) :
        AirTemperature(airTemperature),
        ConvectionCoefficient(convectionCoefficient),
        AirHumidity(airHumidity)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // VariableBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    VariableBCHCCoefficients::VariableBCHCCoefficients(const double airTemperature,
                                                       const double airHumidity) :
        AirTemperature(airTemperature), AirHumidity(airHumidity)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // VariableBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    LinearizedRadiationBCCoefficients::LinearizedRadiationBCCoefficients(
      const double radiationCoefficient, const double radiationTemperature) :
        RadiationCoefficient(radiationCoefficient), RadiationTemperature(radiationTemperature)
    {}
}   // namespace HygroThermFEM
