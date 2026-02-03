#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IEnvironmentCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    IEnvironmentCoefficients::IEnvironmentCoefficients(const double airTemperature,
                                                       const double airHumidity) :
        AirTemperature(airTemperature), AirHumidity(airHumidity)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // FixedBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    FixedBCHCCoefficients::FixedBCHCCoefficients(const double airTemperature,
                                                 const double convectionCoefficient) :
        IEnvironmentCoefficients(airTemperature, 0), ConvectionCoefficient(convectionCoefficient)
    {}

    FixedBCHCCoefficients::FixedBCHCCoefficients(const double airTemperature,
                                                 const double convectionCoefficient,
                                                 const double airHumidity) :
        IEnvironmentCoefficients(airTemperature, airHumidity),
        ConvectionCoefficient(convectionCoefficient)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TemperatureAndHumidity
    ///////////////////////////////////////////////////////////////////////////////////////////////

    TemperatureAndHumidity::TemperatureAndHumidity(const double temperature,
                                                   const double humidity) :
        Temperature(temperature), Humidity(humidity)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TARPCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    TARPCoefficients::TARPCoefficients(const double airTemperature, const double airHumidity) :
        IEnvironmentCoefficients(airTemperature, airHumidity)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// ASHRAEInsideCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ASHRAEInsideCoefficients::ASHRAEInsideCoefficients(const double air_temperature,
                                                       const double air_humidity,
                                                       const double air_pressure) :
        IEnvironmentCoefficients(air_temperature, air_humidity), AirPressure(air_pressure)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// ASHRAEOutsideCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ASHRAEOutsideCoefficients::ASHRAEOutsideCoefficients(const double airTemperature,
                                                         const double airHumidity,
                                                         const double windSpeed) :
        IEnvironmentCoefficients(airTemperature, airHumidity), WindSpeed(windSpeed)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// YazdanianKlemsCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    YazdanianKlemsCoefficients::YazdanianKlemsCoefficients(const double airTemperature,
                                                           const double airHumidity,
                                                           const double windSpeed,
                                                           const WindDirection windDirection) :
        IEnvironmentCoefficients(airTemperature, airHumidity),
        WindSpeed(windSpeed),
        WindDir(windDirection)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// KimuraCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    KimuraCoefficients::KimuraCoefficients(const double airTemperature,
                                           const double airHumidity,
                                           const double windSpeed,
                                           const WindDirection windDir) :
        IEnvironmentCoefficients(airTemperature, airHumidity),
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
