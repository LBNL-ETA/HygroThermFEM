#include "ConvectiveCoefficient.hxx"

#include <algorithm>

#include "Common.hxx"
#include "FEMMath.hxx"
#include "Node2D.hxx"
#include "VectorOperators.hxx"
#include "IBCLine2D.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// IConvectiveCoefficient
    ////////////////////////////////////////////////////////

    IConvectiveCoefficient::IConvectiveCoefficient(const INodes & nodes) : m_Nodes(nodes)
    {}

    std::vector<double> IConvectiveCoefficient::waterVaporTransferCoefficient() const
    {
        std::vector<double> beta(m_Nodes.size(), 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double humidity = m_Nodes[j].property(Variable::humidity);
            if(humidity <= 1)
            {
                beta[j] = 1 / (Constants::Cp_Air * Constants::Density_Air);
            }
            else
            {
                beta[j] = 0;
            }
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

    std::vector<double> FixedConvectionCoefficient::convectiveCoefficients() const
    {
        std::vector<double> result(m_Nodes.size(), m_ConvectionFilmCoefficient);

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

    std::vector<double> TARPFilmCoefficient::convectiveCoefficients() const
    {
        const auto minimumConvectionCoefficient = 3.0;
        std::vector<double> result;
        for(const auto & temperature : m_Nodes.properties(Variable::temperature))
        {
            result.push_back(
              std::max(minimumConvectionCoefficient,
                       1.81 * std::pow(std::abs(temperature - m_AirTemperature), 1.0 / 3.0)
                         / (1.382 - std::abs(std::cos(radians(m_SurfaceTilt))))));
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

    std::vector<double> ASHRAEInsideFilmCoefficient::convectiveCoefficients() const
    {
        std::vector<double> result;
        const auto surfaceTiltRad{radians(m_SurfaceTilt)};
        for(const auto & temperature : m_Nodes.properties(Variable::temperature))
        {
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
            result.push_back(Gnui * (prop.m_ThermalConductivity / m_SurfaceHeight));
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

    std::vector<double> ASHRAEOutsideFilmCoefficient::convectiveCoefficients() const
    {
        std::vector<double> result;
        for(size_t i = 0u; i < m_Nodes.size(); ++i)
        {
            result.push_back(4 + 4 * m_WindSpeed);
        }
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

    std::vector<double> YazdanianKlemsFilmCoefficient::convectiveCoefficients() const
    {
        std::vector<double> result;
        for(const auto & temperature : m_Nodes.properties(Variable::temperature))
        {
            const auto first{
              std::pow(0.84 * std::pow(std::abs(temperature - m_AirTemperature), 0.33), 2)};
            const auto second{std::pow(
              coeffs.at(m_Direction).A * std::pow(m_WindSpeed, coeffs.at(m_Direction).B), 2)};
            result.push_back(std::pow(first + second, 0.5));
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

    std::vector<double> KimuraFilmCoefficient::convectiveCoefficients() const
    {
        std::vector<double> result;
        std::map<WindDirection, double> Vc{
          {WindDirection::Windward, m_WindSpeed > 2 ? 0.25 * m_WindSpeed : 0.5 * m_WindSpeed},
          {WindDirection::Leeward, 0.3 + 0.05 * m_WindSpeed}};
        const auto filmCoefficient{4.7 + 7.6 * Vc.at(m_Direction)};

        for(auto i = 0u; i < m_Nodes.size(); ++i)
        {
            result.push_back(filmCoefficient);
        }
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
