#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_VariableTARPHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state(
      {.temperature = 0.0, .humidity = 0.0, .pressure = 101325.0, .liquidPercent = 1.0});

    TestHelper::SlabBuilder(multiDomain)
      .gridXCoordinates({0, 0.05, 0.1})
      .height(0.05)
      .material(material.name())
      .state(state)
      .startCorner(TestHelper::StartCorner::BottomRight)
      .direction(TestHelper::Direction::CounterClockwise)
      .build();

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::TARPCoefficients> bcCoeff{{20.0, 0.6},
                                                               {20.0, 0.5},
                                                               {20.0, 0.4},
                                                               {20.0, 0.3},
                                                               {20.0, 0.2},
                                                               {18.0, 0.2},
                                                               {16.0, 0.2},
                                                               {14.0, 0.2},
                                                               {12.0, 0.2},
                                                               {10.0, 0.2}};

    const auto surfaceTilt{90.0};

    multiDomain.createBC_TARPHc(1, 2, bcCoeff, surfaceTilt);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{{3.46405340748, 3.46405340748, 0.00336122089849, 0.00336122089849, 6.5228697861e-06, 6.5228697861e-06},
 {5.44163791419, 5.44163791419, 0.01030528727, 0.01030528727, 2.97898834034e-05, 2.97898834034e-05},
 {6.48430357621, 6.48430357621, 0.0188467643956, 0.0188467643956, 7.80861273822e-05, 7.80861273822e-05},
 {6.77713832865, 6.77713832866, 0.028211816818, 0.028211816818, 0.000157610418779, 0.000157610418779},
 {6.45365197927, 6.45365197927, 0.03772917704, 0.03772917704, 0.000271658600155, 0.000271658600155},
 {6.02270303924, 6.02270303924, 0.0471511425172, 0.0471511425172, 0.000421476080728, 0.000421476080728},
 {5.51590147706, 5.51590147704, 0.0562975946496, 0.0562975946496, 0.000606677195267, 0.000606677195267},
 {4.99829405022, 4.9982940502, 0.0648083487852, 0.0648083487852, 0.000824748609174, 0.000824748609174},
 {4.51082793037, 4.51082793038, 0.072478780361, 0.072478780361, 0.00107175513947, 0.00107175513947},
 {4.06069327367, 4.06069327369, 0.0793236650497, 0.0793236650497, 0.00134329891194, 0.00134329891194}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.96229279689, 3.96229279689, 2.06013544715, 2.06013544715, 1.54186822385, 1.54186822385},
 {5.5366244215, 5.5366244215, 3.8013811484, 3.8013811484, 3.23295973116, 3.23295973116},
 {6.49604698721, 6.49604698721, 5.13059908692, 5.13059908692, 4.65322103229, 4.65322103229},
 {7.14838886593, 7.14838886593, 6.12059362452, 6.12059362452, 5.75146798132, 5.75146798132},
 {7.58976293399, 7.58976293399, 6.8399798171, 6.8399798171, 6.5661765689, 6.5661765689},
 {7.88614912026, 7.88614912026, 7.35216501201, 7.35216501201, 7.15448105005, 7.15448105005},
 {8.0572466934, 8.0572466934, 7.69707034728, 7.69707034728, 7.56063209556, 7.56063209556},
 {8.12405455528, 8.12405455528, 7.9052651186, 7.9052651186, 7.81863993421, 7.81863993421},
 {8.09861322601, 8.09861322601, 7.99822646479, 7.99822646479, 7.95313405218, 7.95313405218},
 {7.97964816271, 7.97964816271, 7.9861406483, 7.9861406483, 7.97793403466, 7.97793403466}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
