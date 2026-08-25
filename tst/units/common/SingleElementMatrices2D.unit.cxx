#include <gtest/gtest.h>
#include <stdexcept>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

TEST(TestSingleElementMatrices2D, TestConductionMatrix)
{
    SCOPED_TRACE("Begin Test: Single element isothropic conduction matrix and "
                 "RhoCp matrix.");

    // Create MultiDomain before nodes
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    multiDomain.nodes().createNode({.index = 1, .x = 5, .y = 5});
    multiDomain.nodes().createNode({.index = 2, .x = 5, .y = 0});
    multiDomain.nodes().createNode({.index = 3, .x = 15, .y = 0});
    multiDomain.nodes().createNode({.index = 4, .x = 15, .y = 5});

    auto params = TestHelper::TestMaterial();
    params.density = 1.0;
    params.heatCapacity = 1.0;
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    const HygroThermFEM::ElementThermalLinear2D aElem{multiDomain.nodes(),
                                                      multiDomain.materials(),
                                                      1,
                                                      2,
                                                      3,
                                                      4,
                                                      material.name(),
                                                      multiDomain.physicsOptions()};

    auto condMat = aElem.DDuMatrices();

    std::vector<std::vector<double>> correctCondMat = {
      {0.833333333, -0.583333333, -0.416666667, 0.166666667},
      {-0.583333333, 0.833333333, 0.166666667, -0.416666667},
      {-0.416666667, 0.166666667, 0.833333333, -0.583333333},
      {0.166666667, -0.416666667, -0.583333333, 0.833333333}};

    for(auto i = 0; i < 4; ++i)
    {
        for(auto j = 0; j < 4; ++j)
        {
            EXPECT_NEAR(correctCondMat[i][j], condMat(i, j), 1e-6);
        }
    }

    auto rhoCpMat = aElem.capacitanceMatrices();

    std::vector<std::vector<double>> correctRhoCpMat = {
      {5.555555556, 2.777777778, 1.388888889, 2.777777778},
      {2.777777778, 5.555555556, 2.777777778, 1.388888889},
      {1.388888889, 2.777777778, 5.555555556, 2.777777778},
      {2.777777778, 1.388888889, 2.777777778, 5.555555556}};

    for(auto i = 0; i < 4; ++i)
    {
        for(auto j = 0; j < 4; ++j)
        {
            EXPECT_NEAR(correctRhoCpMat[i][j], rhoCpMat(i, j), 1e-6);
        }
    }
}
