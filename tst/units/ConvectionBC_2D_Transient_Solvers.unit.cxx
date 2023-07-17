#include <memory>
#include <gtest/gtest.h>
#include <fstream>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class ConvectionBC_2D_Transient_Solvers : public testing::Test
{
public:
    HygroThermFEM::SingleDomain domain{HygroThermFEM::DomainType::Thermal};

    const size_t nSteps{20u};
    const double dTime{3600};

    const double hc{1};
    const double temperatureAir{-18.0};

    const double initialT1{10};
    const double initialT2{20};
    const double initialT3{30};

protected:
    void SetUp() override
    {
        HygroThermFEM::State state1;
        state1.setValue(HygroThermFEM::BaseVariable::temperature, initialT1);
        HygroThermFEM::State state2;
        state2.setValue(HygroThermFEM::BaseVariable::temperature, initialT2);
        HygroThermFEM::State state3;
        state3.setValue(HygroThermFEM::BaseVariable::temperature, initialT3);

        NodePool::Instance().createNode(1, 0.2, 0.05, state1);
        NodePool::Instance().createNode(2, 0.2, 0.00, state1);
        NodePool::Instance().createNode(3, 0.1, 0.05, state2);
        NodePool::Instance().createNode(4, 0.1, 0.00, state2);
        NodePool::Instance().createNode(5, 0.0, 0.05, state3);
        NodePool::Instance().createNode(6, 0.0, 0.00, state3);

        // Material Properties
        const double thermalConductivityDry{1.0};
        const double density{2050.0};
        const double porosity{0.0};
        const double specificHeatCapacityDry{850.0};
        const double diffusionResistanceFactor{15.0};
        const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
          {0.0, 1.0}, {180, 1.0}};
        const double thermalConductivityMeasuredAtTemperature{0};
        const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
          {0.0, 1.0}, {1, 1.0}};
        const double thermalConductivityMeasuredAtHumidity{0};
        const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                                  {27, 1E-8},
                                                                                  {45, 1.1E-8},
                                                                                  {90, 2E-8},
                                                                                  {126, 3.5E-8},
                                                                                  {144, 5E-8},
                                                                                  {162, 1E-7},
                                                                                  {171, 2E-7},
                                                                                  {180, 7E-7}};

        const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                                {0.5, 5.3},
                                                                                {0.65, 8.4},
                                                                                {0.8, 12},
                                                                                {0.93, 17},
                                                                                {0.95, 25},
                                                                                {0.99, 63},
                                                                                {0.995, 83},
                                                                                {0.999, 120},
                                                                                {1, 180}};

        auto & material =
          MaterialPool::Instance().createSolidMaterial("Test Material",
                                                       thermalConductivityDry,
                                                       density,
                                                       porosity,
                                                       specificHeatCapacityDry,
                                                       diffusionResistanceFactor,
                                                       thermalConductivityMoistureDependent,
                                                       thermalConductivityMeasuredAtTemperature,
                                                       thermalConductivityTemperatureDependent,
                                                       thermalConductivityMeasuredAtHumidity,
                                                       liquidTransportationCurve,
                                                       moistureStorageFunction);

        createElement(domain, 3, 4, 2, 1, material.name());
        createElement(domain, 6, 4, 3, 5, material.name());

        const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir, hc};
        HygroThermFEM::Thermal::createBC_FixedHc(domain, 6, 5, bcCoeff2);
    }

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(ConvectionBC_2D_Transient_Solvers, Substitution)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    temperaturesSolution.emplace_back(temperatures);

    HygroThermFEM::TransientSubstitutionSolver solver;
    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = solver.transient(domain, temperatures, dTime).solution;
        temperaturesSolution.push_back(temperatures);
    }

    const std::vector<std::vector<double>> correctTemperatureSolution{
      {10.000000, 10.000000, 20.000000, 20.000000, 30.000000, 30.000000},
      {12.864085, 12.864085, 19.795568, 19.795568, 25.737546, 25.737546},
      {14.802374, 14.802374, 19.493304, 19.493304, 22.721188, 22.721188},
      {16.072864, 16.072864, 19.147625, 19.147625, 20.549205, 20.549205},
      {16.866040, 16.866040, 18.785636, 18.785636, 18.953106, 18.953106},
      {17.320493, 17.320493, 18.420331, 18.420331, 17.751993, 17.751993},
      {17.535998, 17.535998, 18.057552, 18.057552, 16.823155, 16.823155},
      {17.583838, 17.583838, 17.699616, 17.699616, 16.082885, 16.082885},
      {17.514642, 17.514642, 17.347179, 17.347179, 15.473819, 15.473819},
      {17.364218, 17.364218, 17.000170, 17.000170, 14.956501, 14.956501},
      {17.157802, 17.157802, 16.658249, 16.658249, 14.503706, 14.503706},
      {16.913138, 16.913138, 16.321018, 16.321018, 14.096603, 14.096603},
      {16.642672, 16.642672, 15.988104, 15.988104, 13.722141, 13.722141},
      {16.355117, 16.355117, 15.659194, 15.659194, 13.371257, 13.371257},
      {16.056567, 16.056567, 15.334036, 15.334036, 13.037649, 13.037649},
      {15.751277, 15.751277, 15.012431, 15.012431, 12.716928, 12.716928},
      {15.442210, 15.442210, 14.694224, 14.694224, 12.406032, 12.406032},
      {15.131429, 15.131429, 14.379296, 14.379296, 12.102823, 12.102823},
      {14.820366, 14.820366, 14.067551, 14.067551, 11.805802, 11.805802},
      {14.510012, 14.510012, 13.758913, 13.758913, 11.513919, 11.513919},
      {14.201050, 14.201050, 13.453319, 13.453319, 11.226433, 11.226433},
    };

    EXPECT_EQ(temperaturesSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperaturesSolution[i][j], 1e-6);
        }
    }
}

TEST_F(ConvectionBC_2D_Transient_Solvers, Sundials)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    temperaturesSolution.emplace_back(temperatures);

    Sundials::SolverIDA solver;
    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution = solver.transient(domain, temperatures, dTime, i);
        temperatures = solution.solution;
        temperaturesSolution.push_back(temperatures);
    }

    const std::vector<std::vector<double>> correctTemperatureSolution{
      {10.000000, 10.000000, 20.000000, 20.000000, 30.000000, 30.000000},
      {13.376549, 13.376549, 19.846799, 19.846799, 25.047837, 25.047837},
      {15.522207, 15.522207, 19.528208, 19.528208, 21.821199, 21.821199},
      {16.807428, 16.807428, 19.160467, 19.160467, 19.680541, 19.680541},
      {17.536859, 17.536859, 18.779016, 18.779016, 18.184592, 18.184592},
      {17.889760, 17.889760, 18.401471, 18.401471, 17.109993, 17.109993},
      {17.996281, 17.996281, 18.028447, 18.028447, 16.298061, 16.298061},
      {17.940890, 17.940890, 17.665353, 17.665353, 15.664700, 15.664700},
      {17.782027, 17.782027, 17.309521, 17.309521, 15.144148, 15.144148},
      {17.567553, 17.567553, 16.974350, 16.974350, 14.714694, 14.714694},
      {17.302581, 17.302581, 16.629831, 16.629831, 14.314515, 14.314515},
      {17.011378, 17.011378, 16.289984, 16.289984, 13.947164, 13.947164},
      {16.704377, 16.704377, 15.954326, 15.954326, 13.602042, 13.602042},
      {16.388322, 16.388322, 15.622515, 15.622515, 13.272314, 13.272314},
      {16.062195, 16.062195, 15.288874, 15.288874, 12.948326, 12.948326},
      {15.747591, 15.747591, 14.972251, 14.972251, 12.645490, 12.645490},
      {15.432879, 15.432879, 14.658787, 14.658787, 12.348615, 12.348615},
      {15.119134, 15.119134, 14.348401, 14.348401, 12.056574, 12.056574},
      {14.783125, 14.783125, 14.017508, 14.017508, 11.746624, 11.746624},
      {14.473316, 14.473316, 13.713324, 13.713324, 11.462526, 11.462526},
      {14.165905, 14.165905, 13.412056, 13.412056, 11.181668, 11.181668},
    };

    EXPECT_EQ(temperaturesSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperaturesSolution[i][j], 1e-6);
        }
    }
}
