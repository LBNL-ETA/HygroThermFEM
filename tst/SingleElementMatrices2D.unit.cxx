#include <gtest/gtest.h>
#include <stdexcept>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class TestSingleElementMatrices2D : public testing::Test {

protected:
  void SetUp() override {}

  void TearDown() override {
    NodePool::Instance().clear();
    MaterialPool::Instance().clear();
  }
};

TEST_F(TestSingleElementMatrices2D, TestConductionMatrix) {
  SCOPED_TRACE("Begin Test: Single element isothropic conduction matrix and "
               "RhoCp matrix.");

  // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

  auto node1 = NodePool::Instance().createNode(1, 5, 5);
  auto node2 = NodePool::Instance().createNode(2, 5, 0);
  auto node3 = NodePool::Instance().createNode(3, 15, 0);
  auto node4 = NodePool::Instance().createNode(4, 15, 5);

  auto &material =
	  MaterialPool::Instance().createMaterial("Test Material",
                                  1,       /// density
                                  0.00,    /// porosity
                                  1,       /// specific heat capacity (dry)
                                  1,       /// thermal conductivity (dry)
                                  15,      /// diffusion resistance factor
                                  {{0, 0}, /// liquid transportation coefficient
                                   {27, 1E-8},
                                   {45, 1.1E-8},
                                   {90, 2E-8},
                                   {126, 3.5E-8},
                                   {144, 5E-8},
                                   {162, 1E-7},
                                   {171, 2E-7},
                                   {180, 7E-7}},
                                  {{0, 0}, /// sorption curve
                                   {0.5, 5.3},
                                   {0.65, 8.4},
                                   {0.8, 12},
                                   {0.93, 17},
                                   {0.95, 25},
                                   {0.99, 63},
                                   {0.995, 83},
                                   {0.999, 120},
                                   {1, 180}});

  const MoisThermFEM::ElementThermalLinear2D aElem{node1, node2, node3, node4, material};

  auto condMat = aElem.DDuMatrices();

  std::vector<std::vector<double>> correctCondMat = {
      {0.833333333, -0.583333333, -0.416666667, 0.166666667},
      {-0.583333333, 0.833333333, 0.166666667, -0.416666667},
      {-0.416666667, 0.166666667, 0.833333333, -0.583333333},
      {0.166666667, -0.416666667, -0.583333333, 0.833333333}};

  for (auto i = 0; i < 4; ++i) {
    for (auto j = 0; j < 4; ++j) {
      EXPECT_NEAR(correctCondMat[i][j], condMat(i, j), 1e-6);
    }
  }

  auto rhoCpMat = aElem.capacitanceMatrices();

  std::vector<std::vector<double>> correctRhoCpMat = {
      {5.555555556, 2.777777778, 1.388888889, 2.777777778},
      {2.777777778, 5.555555556, 2.777777778, 1.388888889},
      {1.388888889, 2.777777778, 5.555555556, 2.777777778},
      {2.777777778, 1.388888889, 2.777777778, 5.555555556}};

  for (auto i = 0; i < 4; ++i) {
    for (auto j = 0; j < 4; ++j) {
      EXPECT_NEAR(correctRhoCpMat[i][j], rhoCpMat(i, j), 1e-6);
    }
  }
}
