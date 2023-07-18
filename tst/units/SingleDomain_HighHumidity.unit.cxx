#include <memory>
#include <gtest/gtest.h>
#include <chrono>
#include <fstream>

#include "HygroThermFEM2D.hxx"

class SingleDomain_HighHumidity : public testing::Test
{
public:
    HygroThermFEM::SingleDomain domain{HygroThermFEM::DomainType::Moisture};

    const double dTime{3600};
    const size_t nSteps{20u};

    const double domainTemperature{25.0};
    const double domainHumidity{0.9};
    const double domainPressure{101325.0};
    const double liquidPercent{1.0};

    const double hc{8.0};
    const double airTemperature{15.0};
    const double airHumidity{1.0};

protected:
    void SetUp() override
    {
        using HygroThermFEM::NodePool;
        using HygroThermFEM::MaterialPool;

        std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

        HygroThermFEM::State state(
          domainTemperature, domainHumidity, domainPressure, liquidPercent);
        size_t nodeIndex = 0;
        for(auto val : gridXCoordinates)
        {
            ++nodeIndex;
            NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
            ++nodeIndex;
            NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        }

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

        /// Create elements
        for(size_t i = 1; i <= (HygroThermFEM::maxNodeIndex() - 2) / 2; ++i)
        {
            const auto node1 = 2u * i + 1u;
            const auto node2 = 2u * i + 2u;
            const auto node3 = 2u * i;
            const auto node4 = 2u * i - 1u;
            createElement(domain, node2, node3, node4, node1, material.name());
        }

        const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

        HygroThermFEM::Moisture::createBC_FixedHc(domain, 5, 6, bcCoeff);
    }

    void TearDown() override
    {
        using HygroThermFEM::NodePool;
        using HygroThermFEM::MaterialPool;

        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(SingleDomain_HighHumidity, Substitution)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer (Substitution).");

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> humiditySolution;
    std::vector<double> humidityError;

    humiditySolution.emplace_back(humidities);

    HygroThermFEM::TransientSingleDomainSubstitution solver{domain};
    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution{solver.transient(humidities, dTime, i)};
        humidities = solution.solution;
        humiditySolution.emplace_back(humidities);
    }

    const std::vector<std::vector<double>> correctHumiditySolution{
      {0.900000, 0.900000, 0.900000, 0.900000, 0.900000, 0.900000},
      {0.900000, 0.900000, 0.899984, 0.899984, 0.815263, 0.815263},
      {0.900000, 0.900000, 0.899956, 0.899956, 0.746663, 0.746663},
      {0.900000, 0.900000, 0.899917, 0.899917, 0.694754, 0.694754},
      {0.900000, 0.900000, 0.899872, 0.899872, 0.656811, 0.656811},
      {0.900000, 0.900000, 0.899821, 0.899821, 0.629250, 0.629250},
      {0.900000, 0.900000, 0.899767, 0.899767, 0.609326, 0.609326},
      {0.900000, 0.900000, 0.899710, 0.899710, 0.594974, 0.594974},
      {0.900000, 0.900000, 0.899650, 0.899650, 0.584664, 0.584664},
      {0.900000, 0.900000, 0.899590, 0.899590, 0.577272, 0.577272},
      {0.900000, 0.900000, 0.899528, 0.899528, 0.571979, 0.571979},
      {0.900000, 0.900000, 0.899466, 0.899466, 0.568193, 0.568193},
      {0.900000, 0.900000, 0.899403, 0.899403, 0.565489, 0.565489},
      {0.899999, 0.899999, 0.899339, 0.899339, 0.563560, 0.563560},
      {0.899999, 0.899999, 0.899275, 0.899275, 0.562183, 0.562183},
      {0.899999, 0.899999, 0.899211, 0.899211, 0.561203, 0.561203},
      {0.899999, 0.899999, 0.899147, 0.899147, 0.560504, 0.560504},
      {0.899999, 0.899999, 0.899083, 0.899083, 0.560007, 0.560007},
      {0.899999, 0.899999, 0.899018, 0.899018, 0.559653, 0.559653},
      {0.899999, 0.899999, 0.898953, 0.898953, 0.559401, 0.559401},
      {0.899999, 0.899999, 0.898889, 0.898889, 0.559221, 0.559221},
    };

    EXPECT_EQ(humiditySolution.size(), correctHumiditySolution.size());

    for(auto i = 0u; i < humiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < humiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][j], 1e-6);
        }
    }
}

TEST_F(SingleDomain_HighHumidity, Sundials)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer (Sundials).");

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> humiditySolution;

    humiditySolution.emplace_back(humidities);

    Sundials::TransientSingleDomainSundials solver{domain};
    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution{solver.transient(humidities, dTime, i)};
        humidities = solution.solution;
        humiditySolution.emplace_back(humidities);
    }

    const std::vector<std::vector<double>> correctHumiditySolution{
      {0.900000, 0.900000, 0.900000, 0.900000, 0.900000, 0.900000},
      {0.900000, 0.900000, 0.899997, 0.899997, 0.866514, 0.866514},
      {0.900000, 0.900000, 0.899983, 0.899983, 0.791222, 0.791222},
      {0.900000, 0.900000, 0.899956, 0.899956, 0.723636, 0.723636},
      {0.900000, 0.900000, 0.899918, 0.899918, 0.671543, 0.671543},
      {0.900000, 0.900000, 0.899873, 0.899873, 0.636146, 0.636146},
      {0.900000, 0.900000, 0.899822, 0.899822, 0.611443, 0.611443},
      {0.900000, 0.900000, 0.899765, 0.899765, 0.593958, 0.593958},
      {0.900000, 0.900000, 0.899705, 0.899705, 0.582241, 0.582241},
      {0.900000, 0.900000, 0.899646, 0.899646, 0.574639, 0.574639},
      {0.900000, 0.900000, 0.899585, 0.899585, 0.569492, 0.569492},
      {0.900000, 0.900000, 0.899519, 0.899519, 0.565803, 0.565803},
      {0.900000, 0.900000, 0.899457, 0.899457, 0.563484, 0.563484},
      {0.900000, 0.900000, 0.899397, 0.899397, 0.562008, 0.562008},
      {0.900000, 0.900000, 0.899332, 0.899332, 0.560916, 0.560916},
      {0.899999, 0.899999, 0.899267, 0.899267, 0.560202, 0.560202},
      {0.899999, 0.899999, 0.899199, 0.899199, 0.559708, 0.559708},
      {0.899999, 0.899999, 0.899140, 0.899140, 0.559421, 0.559421},
      {0.899999, 0.899999, 0.899074, 0.899074, 0.559202, 0.559202},
      {0.899999, 0.899999, 0.899010, 0.899010, 0.559064, 0.559064},
      {0.899999, 0.899999, 0.898943, 0.898943, 0.558966, 0.558966},
    };

    EXPECT_EQ(humiditySolution.size(), correctHumiditySolution.size());

    for(auto i = 0u; i < humiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < humiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][j], 1e-6);
        }
    }
}