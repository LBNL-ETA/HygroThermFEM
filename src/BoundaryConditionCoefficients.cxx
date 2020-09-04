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
    // TARPCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    TARPCoefficients::TARPCoefficients(const double airTemperature,
                                                               const double airHumidity) :
        AirTemperature(airTemperature), AirHumidity(airHumidity)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// ASHRAEInsideCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    ASHRAEInsideCoefficients::ASHRAEInsideCoefficients(double air_temperature,
                                                       double air_humidity,
                                                       double air_pressure):
        AirTemperature(air_temperature),
        AirHumidity(air_humidity),
        AirPressure(air_pressure)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// ASHRAEOutsideCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ASHRAEOutsideCoefficients::ASHRAEOutsideCoefficients(double airTemperature,
                                                         double airHumidity,
                                                         double windSpeed) :
      AirTemperature(airTemperature), AirHumidity(airHumidity), WindSpeed(windSpeed)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// YazdanianKlemsCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    YazdanianKlemsCoefficients::YazdanianKlemsCoefficients(double airTemperature,
                                                           double airHumidity,
                                                           double windSpeed,
                                                           WindDirection windDirection) :
      AirTemperature(airTemperature),
      AirHumidity(airHumidity),
      WindSpeed(windSpeed),
      WindDir(windDirection)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// KimuraCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    KimuraCoefficients::KimuraCoefficients(double airTemperature,
                                           double airHumidity,
                                           double windSpeed,
                                           WindDirection windDir) :
      AirTemperature(airTemperature),
      AirHumidity(airHumidity),
      WindSpeed(windSpeed),
      WindDir(windDir)
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
