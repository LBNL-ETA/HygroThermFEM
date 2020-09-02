#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class Moisture_2D_TwoElements_1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(Moisture_2D_TwoElements_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.9;
    const auto initialPressure = 101325.0;
    const auto liquidPercent = 1.0;

    const State state(initialTemperature, initialHumidity, initialPressure, liquidPercent);
    NodePool::Instance().createNode(1, 0.15, 0.05, state);
    NodePool::Instance().createNode(2, 0.15, 0, state);
    NodePool::Instance().createNode(3, 0.05, 0.05, state);
    NodePool::Instance().createNode(4, 0.05, 0, state);
    NodePool::Instance().createNode(5, 0, 0.05, state);
    NodePool::Instance().createNode(6, 0, 0, state);

    // Material Properties (Cottaer Sandstone)
    const double thermalConductivityDry{1.8};
    const double density{2050.0};
    const double porosity{0.22};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
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
      MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone",
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

    HygroThermFEM::MoistureDomain domain;

    /// Create elements
    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto airTemperature = 20.0;
    const auto airHumidity = 0.0;
    const auto hc = 10.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    domain.createBC_FixedHc(5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 24;

    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<double> timesteps;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution = domain.transient(humidities, dTime);
        humidities = solution.solution;
        timesteps.push_back(solution.dTime);
        auto waterContent = NodePool::Instance().properties(HygroThermFEM::Variable::water);
        waterContentSolution.push_back(waterContent);
        fluxSolution.push_back(domain.flux());
    }

    std::vector<double> correctTimesteps{3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600};

    std::vector<std::vector<double>> correctWaterContentSolution{
      {15.8461537, 15.8461537, 15.8448713, 15.8448713, 8.68767655, 8.68767655},
      {15.8461534, 15.8461534, 15.8425833, 15.8425833, 5.0504068, 5.0504068},
      {15.8461527, 15.8461527, 15.8395261, 15.8395261, 3.55745897, 3.55745897},
      {15.8461517, 15.8461517, 15.8359046, 15.8359046, 2.47108274, 2.47108274},
      {15.8461503, 15.8461503, 15.8318899, 15.8318899, 1.71931163, 1.71931163},
      {15.8461484, 15.8461484, 15.8275991, 15.8275991, 1.19935, 1.19935},
      {15.8461461, 15.8461461, 15.8231134, 15.8231134, 0.839926159, 0.839926159},
      {15.8461433, 15.8461433, 15.8184889, 15.8184889, 0.591629538, 0.591629538},
      {15.8461401, 15.8461401, 15.8137648, 15.8137648, 0.420214977, 0.420214977},
      {15.8461364, 15.8461364, 15.8089678, 15.8089678, 0.301958047, 0.301958047},
      {15.8461321, 15.8461321, 15.8041167, 15.8041167, 0.220431377, 0.220431377},
      {15.8461274, 15.8461274, 15.7992244, 15.7992244, 0.164266904, 0.164266904},
      {15.8461222, 15.8461222, 15.7942996, 15.7942996, 0.125602682, 0.125602682},
      {15.8461165, 15.8461165, 15.7893486, 15.7893486, 0.0990051665, 0.0990051665},
      {15.8461103, 15.8461103, 15.7843754, 15.7843754, 0.080721732, 0.080721732},
      {15.8461036, 15.8461036, 15.779383, 15.779383, 0.0681625154, 0.0681625154},
      {15.8460964, 15.8460964, 15.7743733, 15.7743733, 0.0595414311, 0.0595414311},
      {15.8460887, 15.8460887, 15.7693476, 15.7693476, 0.0536276271, 0.0536276271},
      {15.8460805, 15.8460805, 15.7643067, 15.7643067, 0.0495735176, 0.0495735176},
      {15.8460717, 15.8460717, 15.7592513, 15.7592513, 0.0467958837, 0.0467958837},
      {15.8460625, 15.8460625, 15.7541816, 15.7541816, 0.0448937232, 0.0448937232},
      {15.8460527, 15.8460527, 15.7490978, 15.7490978, 0.0435915333, 0.0435915333},
      {15.8460424, 15.8460424, 15.7440001, 15.7440001, 0.0427001811, 0.0427001811},
      {15.8460316, 15.8460316, 15.7388886, 15.7388886, 0.0420899289, 0.0420899289}};

    std::vector<std::vector<HygroThermFEM::NodeFlux>> correctFluxSolution{{{-9.57765514e-12, 0},
                                                                           {-9.57765514e-12, 0},
                                                                           {-6.8365037e-08, 0},
                                                                           {-6.8365037e-08, 0},
                                                                           {-1.36720496e-07, 0},
                                                                           {-1.36720496e-07, 0}},
                                                                          {{-2.66629774e-11, 0},
                                                                           {-2.66629774e-11, 0},
                                                                           {-1.21651251e-07, 0},
                                                                           {-1.21651251e-07, 0},
                                                                           {-2.43275839e-07, 0},
                                                                           {-2.43275839e-07, 0}},
                                                                          {{-4.94907273e-11, 0},
                                                                           {-4.94907273e-11, 0},
                                                                           {-1.62097309e-07, 0},
                                                                           {-1.62097309e-07, 0},
                                                                           {-3.24145128e-07, 0},
                                                                           {-3.24145128e-07, 0}},
                                                                          {{-7.65308467e-11, 0},
                                                                           {-7.65308467e-11, 0},
                                                                           {-1.91523942e-07, 0},
                                                                           {-1.91523942e-07, 0},
                                                                           {-3.82971352e-07, 0},
                                                                           {-3.82971352e-07, 0}},
                                                                          {{-1.06504118e-10, 0},
                                                                           {-1.06504118e-10, 0},
                                                                           {-2.1188144e-07, 0},
                                                                           {-2.1188144e-07, 0},
                                                                           {-4.23656376e-07, 0},
                                                                           {-4.23656376e-07, 0}},
                                                                          {{-1.38536102e-10, 0},
                                                                           {-1.38536102e-10, 0},
                                                                           {-2.25956037e-07, 0},
                                                                           {-2.25956037e-07, 0},
                                                                           {-4.51773539e-07, 0},
                                                                           {-4.51773539e-07, 0}},
                                                                          {{-1.72020592e-10, 0},
                                                                           {-1.72020592e-10, 0},
                                                                           {-2.35679446e-07, 0},
                                                                           {-2.35679446e-07, 0},
                                                                           {-4.71186871e-07, 0},
                                                                           {-4.71186871e-07, 0}},
                                                                          {{-2.06537477e-10, 0},
                                                                           {-2.06537477e-10, 0},
                                                                           {-2.42390865e-07, 0},
                                                                           {-2.42390865e-07, 0},
                                                                           {-4.84575192e-07, 0},
                                                                           {-4.84575192e-07, 0}},
                                                                          {{-2.41795753e-10, 0},
                                                                           {-2.41795753e-10, 0},
                                                                           {-2.47018456e-07, 0},
                                                                           {-2.47018456e-07, 0},
                                                                           {-4.93795116e-07, 0},
                                                                           {-4.93795116e-07, 0}},
                                                                          {{-2.77594001e-10, 0},
                                                                           {-2.77594001e-10, 0},
                                                                           {-2.50205236e-07, 0},
                                                                           {-2.50205236e-07, 0},
                                                                           {-5.00132877e-07, 0},
                                                                           {-5.00132877e-07, 0}},
                                                                          {{-3.13792944e-10, 0},
                                                                           {-3.13792944e-10, 0},
                                                                           {-2.52396443e-07, 0},
                                                                           {-2.52396443e-07, 0},
                                                                           {-5.04479092e-07, 0},
                                                                           {-5.04479092e-07, 0}},
                                                                          {{-3.50296404e-10, 0},
                                                                           {-3.50296404e-10, 0},
                                                                           {-2.53900192e-07, 0},
                                                                           {-2.53900192e-07, 0},
                                                                           {-5.07450087e-07, 0},
                                                                           {-5.07450087e-07, 0}},
                                                                          {{-3.8703808e-10, 0},
                                                                           {-3.8703808e-10, 0},
                                                                           {-2.5492957e-07, 0},
                                                                           {-2.5492957e-07, 0},
                                                                           {-5.09472102e-07, 0},
                                                                           {-5.09472102e-07, 0}},
                                                                          {{-4.23972372e-10, 0},
                                                                           {-4.23972372e-10, 0},
                                                                           {-2.55631846e-07, 0},
                                                                           {-2.55631846e-07, 0},
                                                                           {-5.1083972e-07, 0},
                                                                           {-5.1083972e-07, 0}},
                                                                          {{-4.61068012e-10, 0},
                                                                           {-4.61068012e-10, 0},
                                                                           {-2.56108729e-07, 0},
                                                                           {-2.56108729e-07, 0},
                                                                           {-5.1175639e-07, 0},
                                                                           {-5.1175639e-07, 0}},
                                                                          {{-4.98303658e-10, 0},
                                                                           {-4.98303658e-10, 0},
                                                                           {-2.56430414e-07, 0},
                                                                           {-2.56430414e-07, 0},
                                                                           {-5.12362524e-07, 0},
                                                                           {-5.12362524e-07, 0}},
                                                                          {{-5.35664826e-10, 0},
                                                                           {-5.35664826e-10, 0},
                                                                           {-2.5664531e-07, 0},
                                                                           {-2.5664531e-07, 0},
                                                                           {-5.12754955e-07, 0},
                                                                           {-5.12754955e-07, 0}},
                                                                          {{-5.73141783e-10, 0},
                                                                           {-5.73141783e-10, 0},
                                                                           {-2.56786778e-07, 0},
                                                                           {-2.56786778e-07, 0},
                                                                           {-5.13000414e-07, 0},
                                                                           {-5.13000414e-07, 0}},
                                                                          {{-6.10728072e-10, 0},
                                                                           {-6.10728072e-10, 0},
                                                                           {-2.5687779e-07, 0},
                                                                           {-2.5687779e-07, 0},
                                                                           {-5.13144851e-07, 0},
                                                                           {-5.13144851e-07, 0}},
                                                                          {{-6.48419503e-10, 0},
                                                                           {-6.48419503e-10, 0},
                                                                           {-2.56934153e-07, 0},
                                                                           {-2.56934153e-07, 0},
                                                                           {-5.13219886e-07, 0},
                                                                           {-5.13219886e-07, 0}},
                                                                          {{-6.86213455e-10, 0},
                                                                           {-6.86213455e-10, 0},
                                                                           {-2.56966735e-07, 0},
                                                                           {-2.56966735e-07, 0},
                                                                           {-5.13247257e-07, 0},
                                                                           {-5.13247257e-07, 0}},
                                                                          {{-7.24108386e-10, 0},
                                                                           {-7.24108386e-10, 0},
                                                                           {-2.56983005e-07, 0},
                                                                           {-2.56983005e-07, 0},
                                                                           {-5.13241901e-07, 0},
                                                                           {-5.13241901e-07, 0}},
                                                                          {{-7.62103506e-10, 0},
                                                                           {-7.62103506e-10, 0},
                                                                           {-2.56988086e-07, 0},
                                                                           {-2.56988086e-07, 0},
                                                                           {-5.13214068e-07, 0},
                                                                           {-5.13214068e-07, 0}},
                                                                          {{-8.00198541e-10, 0},
                                                                           {-8.00198541e-10, 0},
                                                                           {-2.56985496e-07, 0},
                                                                           {-2.56985496e-07, 0},
                                                                           {-5.13170793e-07, 0},
                                                                           {-5.13170793e-07, 0}}};

    EXPECT_EQ(correctTimesteps.size(), timesteps.size());

    for(size_t i = 0u; i < correctTimesteps.size(); ++i)
    {
        EXPECT_NEAR(correctTimesteps[i], timesteps[i], 1e-6);
    }

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
            EXPECT_NEAR(correctFluxSolution[i][j].x, fluxSolution[i][j].x, 1e-12);
            EXPECT_NEAR(correctFluxSolution[i][j].y, fluxSolution[i][j].y, 1e-12);
        }
    }
}
