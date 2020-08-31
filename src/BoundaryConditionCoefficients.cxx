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
    // VariableBCTARPHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    VariableBCTARPHCCoefficients::VariableBCTARPHCCoefficients(const double airTemperature,
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
    // LinearizedRadiationBCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    LinearizedRadiationBCCoefficients::LinearizedRadiationBCCoefficients(
      const double radiationCoefficient, const double radiationTemperature) :
        RadiationCoefficient(radiationCoefficient), RadiationTemperature(radiationTemperature)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // BlackBodyRadiationBCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    BlackBodyRadiationBCCoefficients::BlackBodyRadiationBCCoefficients(const double emissivity,
                                                                       const double temperature) :
        Emissivity(emissivity), Temperature(temperature)
    {}
}   // namespace HygroThermFEM
