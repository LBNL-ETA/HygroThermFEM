#include <memory>
#include <gtest/gtest.h>
#include <fstream>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class ConvectionBC_2D_Transient_Sundials : public testing::Test
{
public:
    const size_t nSteps{200u};
    const double dTime{360};
    const double hc2{1};
    const double temperatureAir2{-18.0};
    const double initialT1{10};
    const double initialT2{20};
    const double initialT3{30};

protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(ConvectionBC_2D_Transient_Sundials, TestExample_Substitution)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

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

    HygroThermFEM::SingleDomain domain{HygroThermFEM::DomainType::Thermal};

    createElement(domain, 3, 4, 2, 1, material.name());
    createElement(domain, 6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    // const auto hc1 = 2.4;
    // const auto temperatureAir1 = 20.0;

    // const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1};

    // const auto hc2 = 15.0;
    // const auto temperatureAir2 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2};

    // HygroThermFEM::Thermal::createBC_FixedHc(domain, 1, 2, bcCoeff1);
    HygroThermFEM::Thermal::createBC_FixedHc(domain, 6, 5, bcCoeff2);

    // const auto dTime = 3600;
    // const auto nSteps = 100;

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    std::vector<double> time;
    temperaturesSolution.emplace_back(temperatures);
    time.emplace_back(0.0);

    for(unsigned i = 0; i < nSteps; ++i)
    {
        HygroThermFEM::TransientSubstitutionSolver solver;
        temperatures = solver.transient(domain, temperatures, dTime).solution;
        temperaturesSolution.push_back(temperatures);
        time.push_back((i + 1) * dTime);
    }

    std::ostringstream filename;
    filename << "Substitution_" << dTime << ".csv";

    std::ofstream myfile;
    myfile.open(filename.str());

    // loop over both the time and temperaturesSolution simultaneously
    for(size_t i = 0; i < temperaturesSolution.size(); i++)
    {
        // first print the time
        myfile << time[i] << ", ";
        // then print the temperatures for this time
        for(auto & val : temperaturesSolution[i])
        {
            myfile << val << ", ";
        }
        myfile << std::endl;
    }
    myfile.close();

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.132130505, 1.132130505, -0.656245111, -0.656245111, -5.621029352, -5.621029352},
      {1.657776414, 1.657776414, -1.47222114, -1.47222114, -8.551769334, -8.551769334},
      {1.788809043, 1.788809043, -2.264759626, -2.264759626, -10.15443472, -10.15443472},
      {1.680631717, 1.680631717, -2.977820826, -2.977820826, -11.08768765, -11.08768765}};

    // EXPECT_EQ(temperaturesSolution.size(), correctTemperatureSolution.size());
    //
    // for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    //{
    //    for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
    //    {
    //        EXPECT_NEAR(correctTemperatureSolution[i][j], temperaturesSolution[i][j], 1e-6);
    //    }
    //}
}

TEST_F(ConvectionBC_2D_Transient_Sundials, TestExample_Sundials)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

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

    HygroThermFEM::SingleDomain domain{HygroThermFEM::DomainType::Thermal};

    createElement(domain, 3, 4, 2, 1, material.name());
    createElement(domain, 6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    // const auto hc1 = 2.4;
    // const auto temperatureAir1 = 20.0;

    // const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1};

    // const auto hc2 = 15.0;
    // const auto temperatureAir2 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2};

    // HygroThermFEM::Thermal::createBC_FixedHc(domain, 1, 2, bcCoeff1);
    HygroThermFEM::Thermal::createBC_FixedHc(domain, 6, 5, bcCoeff2);

    // const auto dTime = 3600;
    // const auto nSteps = 100;

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    std::vector<double> time;
    temperaturesSolution.emplace_back(temperatures);
    time.emplace_back(0.0);


    auto aSolution{Sundials::transient(domain, temperatures, dTime, nSteps)};
    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperaturesSolution.emplace_back(aSolution[i].solution);
        time.emplace_back(aSolution[i].dTime);
    }

    std::ostringstream filename;
    filename << "Sundials_" << dTime << ".csv";

    std::ofstream myfile;
    myfile.open(filename.str());

    // loop over both the time and temperaturesSolution simultaneously
    for(size_t i = 0; i < temperaturesSolution.size(); i++)
    {
        // first print the time
        myfile << time[i] << ", ";
        // then print the temperatures for this time
        for(auto & val : temperaturesSolution[i])
        {
            myfile << val << ", ";
        }
        myfile << std::endl;
    }
    myfile.close();

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.132130505, 1.132130505, -0.656245111, -0.656245111, -5.621029352, -5.621029352},
      {1.657776414, 1.657776414, -1.47222114, -1.47222114, -8.551769334, -8.551769334},
      {1.788809043, 1.788809043, -2.264759626, -2.264759626, -10.15443472, -10.15443472},
      {1.680631717, 1.680631717, -2.977820826, -2.977820826, -11.08768765, -11.08768765}};

    // EXPECT_EQ(temperaturesSolution.size(), correctTemperatureSolution.size());
    //
    // for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    //{
    //    for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
    //    {
    //        EXPECT_NEAR(correctTemperatureSolution[i][j], temperaturesSolution[i][j], 1e-6);
    //    }
    //}
}
