#include <algorithm>

#include "BoundaryCondition2D.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"
#include "VectorOperators.hxx"
#include "Common.hxx"
#include "SimulationProperties.hxx"
#include "FEMMath.hxx"

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
                         / (1.382 - std::abs(std::cos(HygroThermFEM::radians(m_SurfaceTilt))))));
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

    ////////////////////////////////////////////////////////
    /// IConvectionBC
    ////////////////////////////////////////////////////////

    IConvectionBC::IConvectionBC(size_t index1,
                                 size_t index2,
                                 const double t_AirTemperature,
                                 std::unique_ptr<IConvectiveCoefficient> convModel,
                                 const double t_AirHumidity,
                                 const bool t_SimulateMoisture) :
        IBCLinear2D(index1, index2),
        m_AirTemperature(t_AirTemperature),
        m_ConvectiveCoeffCalc(std::move(convModel)),
        m_AirHumidity(t_AirHumidity),
        m_SimulateVaporFluxEnergy(t_SimulateMoisture)
    {}

    std::vector<double> IConvectionBC::R_Vector() const
    {
        auto rightHandSide = m_ConvectiveCoeffCalc->convectiveCoefficients() * m_AirTemperature;

        const auto excludeHeatOfEvaporation =
          SimulationProperties::Instance().excludeHeatOfEvaporation();
        // Moisture is dumping some energy into domain. However, it is possible that user
        // choose not to simulate moisture in which case energy should not be included in
        // simulation
        if(m_SimulateVaporFluxEnergy && !excludeHeatOfEvaporation)
        {
            // Vapor leaking part is added here
            std::vector<double> vaporLeak(numOfBCNodes, 0);
            for(std::size_t j = 0; j < numOfBCNodes; ++j)
            {
                const double T = m_Nodes[j].property(Variable::temperature);
                const double humidity = m_Nodes[j].property(Variable::humidity);
                vaporLeak[j] =
                  (m_AirHumidity * saturationConcentrationAtTemperature(m_AirTemperature)
                   - humidity * saturationConcentrationAtTemperature(T))
                  * heatOfEvaporation(T);
            }
            const auto vaporFluxEnergy =
              vaporLeak * m_ConvectiveCoeffCalc->waterVaporTransferCoefficient();
            rightHandSide = rightHandSide + vaporFluxEnergy;
        }

        return m_PsiVector * rightHandSide;
    }

    SquareMatrix IConvectionBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(m_ConvectiveCoeffCalc->convectiveCoefficients());
    }

    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    ConstantConvectionBC::ConstantConvectionBC(size_t index1,
                                               size_t index2,
                                               const FixedBCHCCoefficients & fixedBCHCCoefficients,
                                               const bool simulateVaporFluxEnergy) :
        IConvectionBC(index1,
                      index2,
                      fixedBCHCCoefficients.AirTemperature,
                      ConvectionModelFactory::createFixedFilmCoefficient(
                        m_Nodes, fixedBCHCCoefficients.ConvectionCoefficient),
                      fixedBCHCCoefficients.AirHumidity,
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// TARPConvectionBC
    ////////////////////////////////////////////////////////
    TARPConvectionBC::TARPConvectionBC(size_t index1,
                                       size_t index2,
                                       const VariableBCTARPHCCoefficients & varHCCoeff,
                                       const double surfaceTilt,
                                       const bool simulateVaporFluxEnergy) :
        IConvectionBC(index1,
                      index2,
                      varHCCoeff.AirTemperature,
                      ConvectionModelFactory::createTARPFilmCoefficient(
                        m_Nodes, varHCCoeff.AirTemperature, surfaceTilt),
                      varHCCoeff.AirHumidity,
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_NodeTemperatures) :
        ConstantConvectionBC(index1, index2, {t_NodeTemperatures, 1e18})
    {
        auto & node1 = NodePool::Instance().getNode(index1);
        auto & node2 = NodePool::Instance().getNode(index2);
        node1.setStateProperty(BaseVariable::temperature, t_NodeTemperatures);
        node2.setStateProperty(BaseVariable::temperature, t_NodeTemperatures);
    }

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_Temp1,
                                 const double t_Temp2) :
        ConstantConvectionBC(index1, index2, {(t_Temp1 + t_Temp2) / 2, 1e18})
    {
        auto & node1 = NodePool::Instance().getNode(index1);
        auto & node2 = NodePool::Instance().getNode(index2);
        node1.setStateProperty(BaseVariable::temperature, t_Temp1);
        node2.setStateProperty(BaseVariable::temperature, t_Temp2);
    }

    ////////////////////////////////////////////////////////
    /// Flux BC
    ////////////////////////////////////////////////////////

    FluxBC::FluxBC(const size_t index1, const size_t index2, const double t_Flux) :
        IBCLinear2D(index1, index2), m_Flux(t_Flux)
    {}

    std::vector<double> FluxBC::R_Vector() const
    {
        std::vector<double> result(m_PsiVector.size(), 0);
        std::transform(m_PsiVector.begin(), m_PsiVector.end(), result.begin(), [&](auto && data) {
            return data * m_Flux;
        });
        return result;
    }

    SquareMatrix FluxBC::H_Matrix() const
    {
        // Flux boundary conditions do not have H matrix (It is zero)
        return SquareMatrix(4);
    }

    ////////////////////////////////////////////////////////
    /// BlackBodyRadiationBC
    ////////////////////////////////////////////////////////

    BlackBodyRadiationBC::BlackBodyRadiationBC(const size_t index1,
                                               const size_t index2,
                                               const double t_Emissivity,
                                               const double t_RadiationTemperature) :
        IBCLinear2D(index1, index2, false),
        m_RadiationTemperature{t_RadiationTemperature},
        m_Emissivity{t_Emissivity}
    {}

    std::vector<double> BlackBodyRadiationBC::HRadiative() const
    {
        std::vector<double> result(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double T = celsiusToKelvin(m_Nodes[j].property(Variable::temperature));
            const double Trad = celsiusToKelvin(m_RadiationTemperature);
            result[j] =
              (T + Trad) * (Trad * Trad + T * T) * Constants::STEFANBOLTZMANN * m_Emissivity;
        }
        return result;
    }

    std::vector<double> BlackBodyRadiationBC::R_Vector() const
    {
        return m_PsiVector * HRadiative() * m_RadiationTemperature;
    }

    SquareMatrix BlackBodyRadiationBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(HRadiative());
    }

    ////////////////////////////////////////////////////////
    /// SimplifiedRadiationBC
    ////////////////////////////////////////////////////////

    LinearizedRadiationBC::LinearizedRadiationBC(
      size_t index1, size_t index2, const LinearizedRadiationBCCoefficients & linearRadBC) :
        IBCLinear2D(index1, index2),
        m_RadiationCoefficient(linearRadBC.RadiationCoefficient),
        m_RadiationTemperature(linearRadBC.RadiationTemperature)
    {}

    std::vector<double> LinearizedRadiationBC::R_Vector() const
    {
        return m_PsiVector * m_RadiationCoefficient * m_RadiationTemperature;
    }

    SquareMatrix LinearizedRadiationBC::H_Matrix() const
    {
        return m_PsiPsiMatrix * m_RadiationCoefficient;
    }

    /////////////////////////////////////////////////////
    /// IMoistureBC
    /////////////////////////////////////////////////////

    IMoistureBC::IMoistureBC(const size_t index1,
                             const size_t index2,
                             const std::string & materialName,
                             const double t_AirHumidity,
                             const double t_AirTemperature,
                             std::unique_ptr<IConvectiveCoefficient> model) :
        IBCLinear2D(index1, index2),
        m_AirHumidity(t_AirHumidity),
        m_AirTemperature(t_AirTemperature),
        m_Material(MaterialPool::Instance().material(materialName)),
        m_ConvectiveCoeffCalc(std::move(model))
    {}

    std::vector<double> IMoistureBC::R_Vector() const
    {
        const auto satOutside = saturationConcentrationAtTemperature(m_AirTemperature);
        const auto gconv =
          m_ConvectiveCoeffCalc->waterVaporTransferCoefficient() * satOutside * m_AirHumidity;
        return m_PsiVector * gconv;
    }

    SquareMatrix IMoistureBC::H_Matrix() const
    {
        std::vector<double> concentration(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double T = m_Nodes[j].property(Variable::temperature);
            concentration[j] = saturationConcentrationAtTemperature(T);
        }
        const auto vaporFlux =
          concentration * m_ConvectiveCoeffCalc->waterVaporTransferCoefficient();

        return m_PsiPsiMatrix.mmultRows(vaporFlux);
    }

    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////
    MoistureBCTARPHc::MoistureBCTARPHc(size_t index1,
                                       size_t index2,
                                       const std::string & materialName,
                                       const VariableBCTARPHCCoefficients & varHCCoeff,
                                       double surfaceTilt) :
        IMoistureBC(index1,
                    index2,
                    materialName,
                    varHCCoeff.AirHumidity,
                    varHCCoeff.AirTemperature,
                    ConvectionModelFactory::createTARPFilmCoefficient(
                      m_Nodes, varHCCoeff.AirTemperature, surfaceTilt))
    {}


    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////
    MoistureBCFixedHc::MoistureBCFixedHc(size_t index1,
                                         size_t index2,
                                         const std::string & materialName,
                                         const FixedBCHCCoefficients & fixedBchcCoefficients) :
        IMoistureBC(index1,
                    index2,
                    materialName,
                    fixedBchcCoefficients.AirHumidity,
                    fixedBchcCoefficients.AirTemperature,
                    ConvectionModelFactory::createFixedFilmCoefficient(
                      m_Nodes, fixedBchcCoefficients.ConvectionCoefficient))
    {}

}   // namespace HygroThermFEM
