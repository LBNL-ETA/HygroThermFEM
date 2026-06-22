#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestVectorOperators, TestAddition)
{
    SCOPED_TRACE("Begin Test: Test Vector Operator - Addition.");

    using HygroThermFEM::operator+;

    const std::vector<double> v1{1, 2, 3};
    const std::vector<double> v2{4, 9, 3};

    const std::vector<double> correctSolution{5, 11, 6};

    const auto solution = v1 + v2;

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        EXPECT_NEAR(correctSolution[i], solution[i], 1e-6);
    }
}

TEST(TestVectorOperators, TestSubtration)
{
    SCOPED_TRACE("Begin Test: Test Vector Operator - Subtraction.");

    using HygroThermFEM::operator-;

    const std::vector<double> v1{1, 2, 3};
    const std::vector<double> v2{4, 9, 3};

    const std::vector<double> correctSolution1{-3, -7, 0};

    const auto solution1 = v1 - v2;

    for(auto i = 0u; i < correctSolution1.size(); ++i)
    {
        EXPECT_NEAR(correctSolution1[i], solution1[i], 1e-6);
    }

    const std::vector<double> correctSolution2{3, 7, 0};

    const auto solution2 = v2 - v1;

    for(auto i = 0u; i < correctSolution2.size(); ++i)
    {
        EXPECT_NEAR(correctSolution2[i], solution2[i], 1e-6);
    }
}

TEST(TestVectorOperators, TestScalarAddition)
{
    SCOPED_TRACE("Begin Test: Test Vector Operator - Scalar Addition.");

    using HygroThermFEM::operator+;

    const std::vector<double> v1{1, 2, 3};
    const auto value = 4.0;

    const std::vector<double> correctSolution{5, 6, 7};

    const auto solution1 = v1 + value;
    const auto solution2 = value + v1;

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        EXPECT_NEAR(correctSolution[i], solution1[i], 1e-6);
        EXPECT_NEAR(correctSolution[i], solution2[i], 1e-6);
    }
}

TEST(TestVectorOperators, TestScalarSubtraction)
{
    SCOPED_TRACE("Begin Test: Test Vector Operator - Scalar Subtraction.");

    using HygroThermFEM::operator-;

    const std::vector<double> v1{1, 2, 3};
    const auto value = 4.0;

    const std::vector<double> correctVectorMinusScalar{-3, -2, -1};
    const auto solution1 = v1 - value;
    for(auto i = 0u; i < correctVectorMinusScalar.size(); ++i)
    {
        EXPECT_NEAR(correctVectorMinusScalar[i], solution1[i], 1e-6);
    }

    const std::vector<double> correctScalarMinusVector{3, 2, 1};
    const auto solution2 = value - v1;
    for(auto i = 0u; i < correctScalarMinusVector.size(); ++i)
    {
        EXPECT_NEAR(correctScalarMinusVector[i], solution2[i], 1e-6);
    }
}

TEST(TestVectorOperators, TestMultiplication)
{
    SCOPED_TRACE("Begin Test: Test Vector Operator - Multiplication.");

    using HygroThermFEM::operator*;

    const std::vector<double> v1{1, 2, 3};
    const std::vector<double> v2{4, 9, 3};

    const std::vector<double> correctSolution1{4, 18, 9};

    const auto solution1 = v1 * v2;

    for(auto i = 0u; i < correctSolution1.size(); ++i)
    {
        EXPECT_NEAR(correctSolution1[i], solution1[i], 1e-6);
    }

    auto const value = 4.0;

    const std::vector<double> correctSolution2{4, 8, 12};

    const auto solution2 = v1 * value;

    for(auto i = 0u; i < correctSolution2.size(); ++i)
    {
        EXPECT_NEAR(correctSolution2[i], solution2[i], 1e-6);
    }

    const std::vector<double> correctSolution3{4, 8, 12};

    const auto solution3 = value * v1;

    for(auto i = 0u; i < correctSolution3.size(); ++i)
    {
        EXPECT_NEAR(correctSolution3[i], solution3[i], 1e-6);
    }
}

TEST(TestVectorOperators, TestDivision)
{
    SCOPED_TRACE("Begin Test: Test Vector Operator - Division.");

    using HygroThermFEM::operator/;

    const std::vector<double> v1{1, 2, 3};
    const std::vector<double> v2{4, 9, 3};

    const std::vector<double> correctSolution1{1.0 / 4.0, 2.0 / 9.0, 1.0};

    const auto solution1 = v1 / v2;

    for(auto i = 0u; i < correctSolution1.size(); ++i)
    {
        EXPECT_NEAR(correctSolution1[i], solution1[i], 1e-6);
    }

    const std::vector<double> correctSolution2{4.0, 9.0 / 2.0, 1.0};

    const auto solution2 = v2 / v1;

    for(auto i = 0u; i < correctSolution2.size(); ++i)
    {
        EXPECT_NEAR(correctSolution2[i], solution2[i], 1e-6);
    }

    auto const value = 4.0;

    const std::vector<double> correctSolution3{1 / 4.0, 2.0 / 4.0, 3.0 / 4.0};

    const auto solution3 = v1 / value;

    for(auto i = 0u; i < correctSolution3.size(); ++i)
    {
        EXPECT_NEAR(correctSolution3[i], solution3[i], 1e-6);
    }

    const std::vector<double> correctSolution4{4.0, 4.0 / 2.0, 4.0 / 3.0};

    const auto solution4 = value / v1;

    for(auto i = 0u; i < correctSolution4.size(); ++i)
    {
        EXPECT_NEAR(correctSolution4[i], solution4[i], 1e-6);
    }
}