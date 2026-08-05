#include "ConvectiveCoefficient.hxx"

#include <algorithm>

#include "Common.hxx"
#include "FEMMath.hxx"
#include "Node2D.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// IConvectiveCoefficient
    ////////////////////////////////////////////////////////

    IConvectiveCoefficient::IConvectiveCoefficient(const INodes & nodes) : m_Nodes(nodes)
    {}

    BCVector IConvectiveCoefficient::waterVaporTransferCoefficient() const
    {
        const double betaValue = 1.0 / (Constants::Cp_Air * Constants::Density_Air);

        // Smooth transition width for tapering beta near saturation (RH = 1.0).
        // Without smoothing, beta drops from full value to zero in a hard step at RH = 1.0,
        // which creates a discontinuity in the Jacobian. During Newton-Raphson iteration the
        // solution can oscillate across that discontinuity: one iterate pushes RH slightly
        // above 1.0 (beta = 0, no vapor transfer), the next iterate pulls it back below 1.0
        // (beta restored, vapor transfer resumes), and so on.
        // A linear ramp over the interval [1.0, 1.0 + transitionWidth] eliminates the
        // discontinuity and allows the solver to converge smoothly near saturation.
        constexpr double transitionWidth = 0.01;

        BCVector beta{};
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double humidity = m_Nodes[j].property(Variable::humidity);

            // smoothFactor = 1.0 when humidity <= 1.0 (full vapor transfer),
            // linearly ramps to 0.0 at humidity = 1.0 + transitionWidth,
            // and stays 0.0 beyond (supersaturation — no vapor transfer).
            const double smoothFactor =
                std::max(0.0, std::min(1.0, (1.0 + transitionWidth - humidity) / transitionWidth));
            beta[j] = betaValue * smoothFactor;
        }

        return convectiveCoefficients() * beta;
    }

    ////////////////////////////////////////////////////////
    /// FixedConvectionCoefficient
    ////////////////////////////////////////////////////////

    FixedConvectionCoefficient::FixedConvectionCoefficient(const INodes & nodes,
                                                           double filmCoefficient) :
        IConvectiveCoefficient(nodes), m_ConvectionFilmCoefficient(filmCoefficient)
    {}

    BCVector FixedConvectionCoefficient::convectiveCoefficients() const
    {
        BCVector result{};
        result.fill(m_ConvectionFilmCoefficient);

        return result;
    }

    ////////////////////////////////////////////////////////
    /// TARPFilmCoefficient
    ////////////////////////////////////////////////////////

    TARPFilmCoefficient::TARPFilmCoefficient(const INodes & nodes,
                                             double airTemperature,
                                             double surfaceTilt) :
        IConvectiveCoefficient(nodes), m_AirTemperature(airTemperature), m_SurfaceTilt(surfaceTilt)
    {}

    BCVector TARPFilmCoefficient::convectiveCoefficients() const
    {
        constexpr auto minimumConvectionCoefficient = 3.0;
        BCVector result{};
        for(std::size_t node = 0; node < numOfBCNodes; ++node)
        {
            const double temperature = m_Nodes[node].property(Variable::temperature);
            result[node] =
              std::max(minimumConvectionCoefficient,
                       1.81 * std::pow(std::abs(temperature - m_AirTemperature), 1.0 / 3.0)
                         / (1.382 - std::abs(std::cos(radians(m_SurfaceTilt)))));
        }
        return result;
    }

    ////////////////////////////////////////////////////////
    /// ASHRAEInsideFilmCoefficient
    ////////////////////////////////////////////////////////

    ASHRAEInsideFilmCoefficient::ASHRAEInsideFilmCoefficient(const INodes & nodes,
                                                             double airTemperature,
                                                             double airPressure,
                                                             double surfaceTilt,
                                                             double mSurfaceHeight) :
        IConvectiveCoefficient(nodes),
        m_AirTemperature(celsiusToKelvin(airTemperature)),   // Note that algorithm require Kelvins
        m_AirPressure(airPressure),
        m_SurfaceTilt(surfaceTilt),
        m_SurfaceHeight(mSurfaceHeight)
    {}

    BCVector ASHRAEInsideFilmCoefficient::convectiveCoefficients() const
    {
        BCVector result{};
        const auto surfaceTiltRad{radians(m_SurfaceTilt)};
        for(std::size_t node = 0; node < numOfBCNodes; ++node)
        {
            const double temperature = m_Nodes[node].property(Variable::temperature);
            const double tMean{m_AirTemperature
                               + 0.25 * (celsiusToKelvin(temperature) - m_AirTemperature)};
            const double deltaT{std::abs(m_AirTemperature - celsiusToKelvin(temperature))};
            Gases::CGas gas;
            gas.setTemperatureAndPressure(tMean, m_AirPressure);
            const auto prop{gas.getGasProperties()};
            const auto gr{Constants::GravityConstant * std::pow(m_SurfaceHeight, 3) * deltaT
                          * std::pow(prop.m_Density, 2) / (tMean * std::pow(prop.m_Viscosity, 2))};
            const auto RaCrit{
              2.5e5 * std::pow((std::exp(0.72 * m_SurfaceTilt) / std::sin(surfaceTiltRad)), 0.2)};
            const auto RaL{gr * prop.m_PrandlNumber};

            auto Gnui{0.0};
            if(m_SurfaceTilt >= 0 && m_SurfaceTilt < 15.0)
            {
                Gnui = 0.13 * std::pow(RaL, 1.0 / 3.0);
            }
            else if(m_SurfaceTilt >= 0 && m_SurfaceTilt <= 90.0)
            {
                if(RaL <= RaCrit)
                {
                    Gnui = 0.56 * std::pow((RaL * std::sin(surfaceTiltRad)), 0.25);
                }
                else
                {
                    Gnui = 0.13 * (std::pow(RaL, 1.0 / 3.0) - std::pow(RaCrit, 1.0 / 3.0))
                           + 0.56 * std::pow((RaCrit * std::sin(surfaceTiltRad)), 0.25);
                }
            }
            else if(m_SurfaceTilt > 90 && m_SurfaceTilt <= 179.0)
            {
                Gnui = 0.56 * std::pow((RaL * std::sin(surfaceTiltRad)), 0.25);
            }
            else if(m_SurfaceTilt > 179 && m_SurfaceTilt <= 180.0)
            {
                Gnui = 0.58 * std::pow(RaL, 1.0 / 3.0);
            }
            result[node] = Gnui * (prop.m_ThermalConductivity / m_SurfaceHeight);
        }
        return result;
    }

    ////////////////////////////////////////////////////////
    /// ASHRAEOutsideFilmCoefficient
    ////////////////////////////////////////////////////////

    ASHRAEOutsideFilmCoefficient::ASHRAEOutsideFilmCoefficient(const INodes & nodes,
                                                               const double windSpeed) :
        IConvectiveCoefficient(nodes), m_WindSpeed(windSpeed)
    {}

    BCVector ASHRAEOutsideFilmCoefficient::convectiveCoefficients() const
    {
        BCVector result{};
        result.fill(4 + 4 * m_WindSpeed);
        return result;
    }

    ////////////////////////////////////////////////////////
    /// YazdanianKlemsFilmCoefficient
    ////////////////////////////////////////////////////////

    YazdanianKlemsFilmCoefficient::YazdanianKlemsFilmCoefficient(const INodes & nodes,
                                                                 double airTemperature,
                                                                 double windSpeed,
                                                                 WindDirection direction) :
        IConvectiveCoefficient(nodes),
        m_AirTemperature(airTemperature),
        m_WindSpeed(windSpeed),
        m_Direction(direction)
    {}

    BCVector YazdanianKlemsFilmCoefficient::convectiveCoefficients() const
    {
        BCVector result{};
        const auto & [coeffA, coeffB] = coeffs.at(m_Direction);
        const auto windPart{std::pow(coeffA * std::pow(m_WindSpeed, coeffB), 2)};
        for(std::size_t node = 0; node < numOfBCNodes; ++node)
        {
            const double temperature = m_Nodes[node].property(Variable::temperature);
            const auto naturalPart{
              std::pow(0.84 * std::pow(std::abs(temperature - m_AirTemperature), 0.33), 2)};
            result[node] = std::pow(naturalPart + windPart, 0.5);
        }
        return result;
    }

    ////////////////////////////////////////////////////////
    /// KimuraFilmCoefficient
    ////////////////////////////////////////////////////////

    KimuraFilmCoefficient::KimuraFilmCoefficient(const INodes & nodes,
                                                 double mWindSpeed,
                                                 WindDirection mDirection) :
        IConvectiveCoefficient(nodes), m_WindSpeed(mWindSpeed), m_Direction(mDirection)
    {}

    BCVector KimuraFilmCoefficient::convectiveCoefficients() const
    {
        const double windwardSpeed{m_WindSpeed > 2 ? 0.25 * m_WindSpeed : 0.5 * m_WindSpeed};
        const double leewardSpeed{0.3 + 0.05 * m_WindSpeed};
        const double correctedSpeed{m_Direction == WindDirection::Windward ? windwardSpeed
                                                                           : leewardSpeed};
        const auto filmCoefficient{4.7 + 7.6 * correctedSpeed};

        BCVector result{};
        result.fill(filmCoefficient);
        return result;
    }

    ////////////////////////////////////////////////////////
    /// ConvectionModelFactory
    ////////////////////////////////////////////////////////

    std::unique_ptr<IConvectiveCoefficient>
      ConvectionModelFactory::createFixedFilmCoefficient(const INodes & nodes,
                                                         double filmCoefficient)
    {
        return std::make_unique<FixedConvectionCoefficient>(nodes, filmCoefficient);
    }

    std::unique_ptr<IConvectiveCoefficient> ConvectionModelFactory::createTARPFilmCoefficient(
      const INodes & nodes, double airTemperature, double surfaceTilt)
    {
        return std::make_unique<TARPFilmCoefficient>(nodes, airTemperature, surfaceTilt);
    }

    std::unique_ptr<IConvectiveCoefficient>
      ConvectionModelFactory::createASHRAEOutsideFilmCoefficient(const INodes & nodes,
                                                                 double windSpeed)
    {
        return std::make_unique<ASHRAEOutsideFilmCoefficient>(nodes, windSpeed);
    }

    std::unique_ptr<IConvectiveCoefficient>
      ConvectionModelFactory::createYazdanianKlemsFilmCoefficient(const INodes & nodes,
                                                                  double airTemperature,
                                                                  double windSpeed,
                                                                  WindDirection direction)
    {
        return std::make_unique<YazdanianKlemsFilmCoefficient>(
          nodes, airTemperature, windSpeed, direction);
    }

    std::unique_ptr<IConvectiveCoefficient> ConvectionModelFactory::createKimuraFilmCoefficient(
      const INodes & nodes, double windSpeed, WindDirection direction)
    {
        return std::make_unique<KimuraFilmCoefficient>(nodes, windSpeed, direction);
    }

    std::unique_ptr<IConvectiveCoefficient>
      ConvectionModelFactory::createASHRAEInsideFilmCoefficient(const INodes & nodes,
                                                                double airTemperature,
                                                                double surfaceTilt,
                                                                double surfaceHeight,
                                                                double airPressure)
    {
        return std::make_unique<ASHRAEInsideFilmCoefficient>(
          nodes, airTemperature, airPressure, surfaceTilt, surfaceHeight);
    }
}
