#pragma once

#include "HygroThermFEM2D.hxx"

namespace TestHelper
{
    //! Common test material parameters
    //! Usage: multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    inline HygroThermFEM::SolidMaterialParams CottaerSandstone()
    {
        return {
            .name = "Cottaer Sandstone",
            .thermalConductivityDry = 1.8,
            .density = 2050.0,
            .porosity = 0.22,
            .heatCapacity = 850.0,
            .diffusionResistanceFactor = 15.0,
            .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
            .temperatureDependentMeasurementHumidity = 0,
            .liquidTransportCurve = {{0, 0},
                                     {27, 1E-8},
                                     {45, 1.1E-8},
                                     {90, 2E-8},
                                     {126, 3.5E-8},
                                     {144, 5E-8},
                                     {162, 1E-7},
                                     {171, 2E-7},
                                     {180, 7E-7}},
            .sorptionCurve = {{0, 0},
                              {0.5, 5.3},
                              {0.65, 8.4},
                              {0.8, 12},
                              {0.93, 17},
                              {0.95, 25},
                              {0.99, 63},
                              {0.995, 83},
                              {0.999, 120},
                              {1, 180}}
        };
    }

    inline HygroThermFEM::SolidMaterialParams CottaerSandstoneNonPorous()
    {
        return {
            .name = "Cottaer Sandstone - non porous",
            .thermalConductivityDry = 1.8,
            .density = 2050.0,
            .porosity = 0.0,
            .heatCapacity = 850.0,
            .diffusionResistanceFactor = 15.0,
            .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
            .temperatureDependentMeasurementHumidity = 0,
            .liquidTransportCurve = {{0, 0},
                                     {27, 1E-8},
                                     {45, 1.1E-8},
                                     {90, 2E-8},
                                     {126, 3.5E-8},
                                     {144, 5E-8},
                                     {162, 1E-7},
                                     {171, 2E-7},
                                     {180, 7E-7}},
            .sorptionCurve = {{0, 0},
                              {0.5, 5.3},
                              {0.65, 8.4},
                              {0.8, 12},
                              {0.93, 17},
                              {0.95, 25},
                              {0.99, 63},
                              {0.995, 83},
                              {0.999, 120},
                              {1, 180}}
        };
    }

    inline HygroThermFEM::SolidMaterialParams ConcreteWC05()
    {
        return {
            .name = "Concrete, w/c=0.5",
            .thermalConductivityDry = 1.6,
            .density = 2300,
            .porosity = 0.18,
            .heatCapacity = 850.0,
            .diffusionResistanceFactor = 92,
            .thermalConductivityMoistureDependent = {{0.0, 1.6}, {150, 1.6}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{0.0, 1.6}, {1, 1.6}},
            .temperatureDependentMeasurementHumidity = 0,
            .liquidTransportCurve = {{0, 0},
                                     {72, 7.4E-11},
                                     {85, 2.5E-10},
                                     {100, 1E-9},
                                     {118, 1.2E-9},
                                     {150, 1.2e-9}},
            .sorptionCurve = {{0, 0},
                              {0.05, 27},
                              {0.1, 32},
                              {0.15, 34},
                              {0.2, 35},
                              {0.3, 37},
                              {0.4, 40},
                              {0.5, 48},
                              {0.6, 58},
                              {0.7, 72},
                              {0.8, 85},
                              {0.9, 100},
                              {0.95, 118},
                              {1, 150}}
        };
    }

    //! Simple test material with thermalConductivityDry = 1.0 and no porosity
    inline HygroThermFEM::SolidMaterialParams TestMaterial()
    {
        return {
            .name = "Test Material",
            .thermalConductivityDry = 1.0,
            .density = 2050.0,
            .porosity = 0.0,
            .heatCapacity = 850.0,
            .diffusionResistanceFactor = 15.0,
            .thermalConductivityMoistureDependent = {{0.0, 1.0}, {180, 1.0}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{0.0, 1.0}, {1, 1.0}},
            .temperatureDependentMeasurementHumidity = 0,
            .liquidTransportCurve = {{0, 0},
                                     {27, 1E-8},
                                     {45, 1.1E-8},
                                     {90, 2E-8},
                                     {126, 3.5E-8},
                                     {144, 5E-8},
                                     {162, 1E-7},
                                     {171, 2E-7},
                                     {180, 7E-7}},
            .sorptionCurve = {{0, 0},
                              {0.5, 5.3},
                              {0.65, 8.4},
                              {0.8, 12},
                              {0.93, 17},
                              {0.95, 25},
                              {0.99, 63},
                              {0.995, 83},
                              {0.999, 120},
                              {1, 180}}
        };
    }

}   // namespace TestHelper
