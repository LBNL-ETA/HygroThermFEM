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
    // ConstantBCTemperatures
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ConstantBCTemperatures::ConstantBCTemperatures(const double temperature1,
                                                   const double temperature2) :
        Temperature1(temperature1), Temperature2(temperature2)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // VariableBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    LinearizedRadiationBCCoefficients::LinearizedRadiationBCCoefficients(
      const double radiationCoefficient, const double radiationTemperature) :
        RadiationCoefficient(radiationCoefficient), RadiationTemperature(radiationTemperature)
    {}
}   // namespace HygroThermFEM
