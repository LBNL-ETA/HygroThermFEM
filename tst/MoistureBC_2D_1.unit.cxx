#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class MoistureBC_2D_1 : public testing::Test
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

TEST_F(MoistureBC_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    std::vector<double> gridXCoordinates{0,
                                         0.015,
                                         0.025,
                                         0.035,
                                         0.045,
                                         0.055,
                                         0.065,
                                         0.075,
                                         0.085,
                                         0.095,
                                         0.105,
                                         0.115,
                                         0.125,
                                         0.135,
                                         0.15};

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.0;
    const auto initialPressure = 101325.0;

    State state(initialTemperature, initialHumidity, initialPressure, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    auto & material = MaterialPool::Instance().createMaterial(
      "Cottaer Sandstone",
      2050,                       /// density
      0.22,                       /// porosity
      850,                        /// specific heat capacity (dry)
      15,                         /// diffusion resistance factor (this is mi value)
      {{0.0, 1.8}, {5.3, 1.8}},   /// thermal conductivity as function of water content
      {{0, 0},                    /// liquid transportation coefficient
                 //{ 27,  1E-8 },
                 //{ 45,  1.1E-8 },
                 //{ 90,  2E-8 },
                 //{ 126, 3.5E-8 },
                 //{ 144, 5E-8 },
                 //{ 162, 1E-7 },
                 //{ 171, 2E-7 },
       {5.3, 7E-7}},
      {{0, 0},   /// sorption curve
                 // { 0.5,   5.3 },
                 // { 0.65,  8.4 },
                 // { 0.8,   12 },
                 // { 0.93,  17 },
                 // { 0.95,  25 },
                 // { 0.99,  63 },
                 // { 0.995, 83 },
                 // { 0.999, 120 },
       {1, 5.3}});

    HygroThermFEM::MoistureDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (HygroThermFEM::NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    // Create Boundary Conditions
    const auto airTemperature = 20;
    const auto airHumidity = 0.5;

    domain.createMoistureBCVariableHc(1, 2, airHumidity, airTemperature);

    const auto dTime = 36000;
    const auto nSteps = 4;
    /// const auto dTime = 60;
    /// const auto nSteps = 3000;

    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> humiditySolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime).solution;
        auto waterContent = NodePool::Instance().properties(HygroThermFEM::Variable::water);
        humiditySolution.push_back(waterContent);
        fluxSolution.push_back(domain.flux());
    }

    // clang-format off
	std::vector<std::vector<double>> correctHumiditySolution{
		{2.5303752,  0.931256186, 0.46178469, 0.228986516, 0.113548229, 0.0563055409, 0.0279204994, 0.013845238,  0.00686592101, 0.00340551739, 0.00169050536, 0.000841909055, 0.000424806945, 0.000225426368, 0.000142984086},
		{2.6003812,  1.37228905,  0.83610797, 0.491774601, 0.282125059, 0.158874306,  0.0881919709, 0.0483998962, 0.0263177325,  0.0142049873,  0.00762717306, 0.00409201248,  0.00222259093,  0.00127456741,  0.000860728015},
		{2.61066881, 1.61113671,  1.09779936, 0.718583764, 0.455612121, 0.281555853,  0.170376182,  0.101317461,  0.059379986,   0.0343875362,  0.0197390215,  0.0112980566,   0.00655032456,  0.00402063722,  0.00286500433},
		{2.61634629, 1.75929546,  1.28284593, 0.901236225, 0.61323934,  0.406029363,  0.262614396,  0.166473245,  0.103725657,   0.0637060574,  0.0387127503,  0.0234438353,   0.0143998598,   0.00937891852,  0.00699666859}};

	std::vector<std::vector<HygroThermFEM::NodeFlux>> correctFluxSolution{
		{{5.77799618e-007, 0.0}, {4.16123197e-007, 0.0}, {1.90310013e-007, 0.0}, {9.43695683e-008, 0.0}, {4.67952982e-008, 0.0}, {2.32044968e-008, 0.0}, {1.15064357e-008, 0.0}, {5.70563883e-009, 0.0}, {2.82908897e-009, 0.0}, {1.4025003e-009, 0.0}, {6.94719362e-010, 0.0}, {3.42995138e-010, 0.0}, {1.6706236e-010, 0.0}, {6.89248734e-011, 0.0}, {2.97883513e-011, 0.0}},
		{{4.43738816e-007, 0.0}, {3.67170611e-007, 0.0}, {2.38613061e-007, 0.0}, {1.50125372e-007, 0.0}, {9.0213578e-008, 0.0}, {5.25544675e-008, 0.0}, {2.9937768e-008, 0.0}, {1.67674722e-008, 0.0}, {9.2665736e-009, 0.0}, {5.06500675e-009, 0.0}, {2.7405432e-009, 0.0}, {1.46460276e-009, 0.0}, {7.63507288e-010, 0.0}, {3.31672526e-010, 0.0}, {1.49529985e-010, 0.0}},
		{{3.61154649e-007, 0.0}, {3.19688042e-007, 0.0}, {2.41875407e-007, 0.0}, {1.74028108e-007, 0.0}, {1.1843141e-007, 0.0}, {7.72968814e-008, 0.0}, {4.88433038e-008, 0.0}, {3.00791682e-008, 0.0}, {1.81375268e-008, 0.0}, {1.07424153e-008, 0.0}, {6.25708236e-009, 0.0}, {3.57404172e-009, 0.0}, {1.97212815e-009, 0.0}, {8.94305798e-010, 0.0}, {4.17557566e-010, 0.0}},
		{{3.09672787e-007, 0.0}, {2.83950781e-007, 0.0}, {2.32527862e-007, 0.0}, {1.81458553e-007, 0.0}, {1.34197487e-007, 0.0}, {9.50168299e-008, 0.0}, {6.49179794e-008, 0.0}, {4.30577019e-008, 0.0}, {2.78491664e-008, 0.0}, {1.76180286e-008, 0.0}, {1.09107717e-008, 0.0}, {6.58861789e-009, 0.0}, {3.81149096e-009, 0.0}, {1.79102056e-009, 0.0}, {8.60763387e-010, 0.0}}
	};
    // clang-format on

    EXPECT_EQ(humiditySolution.size(), correctHumiditySolution.size());

    for(auto i = 0u; i < correctHumiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < correctHumiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][2 * j], 1e-6);
            EXPECT_NEAR(correctFluxSolution[i][j].x, fluxSolution[i][2 * j].x, 1e-12);
            EXPECT_NEAR(correctFluxSolution[i][j].y, fluxSolution[i][2 * j].y, 1e-12);
        }
    }
}
