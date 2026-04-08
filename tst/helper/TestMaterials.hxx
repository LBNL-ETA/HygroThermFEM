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

    // ---------------------------------------------------------------------
    // Materials used in the THERM beam (wall) example.
    // Imported from temp/Materials.xml.
    //
    // Notes on the mapping from the THERM XML to SolidMaterialParams:
    //   - The XML stores two liquid transport curves (Suction and
    //     Redistribution); the codebase uses a single liquidTransportCurve.
    //     Here we pick the curve that has data (Suction first, then
    //     Redistribution).
    //   - thermalConductivityMoistureDependent / temperatureDependent are
    //     duplicated to two points when the XML has only one, so the curves
    //     have a defined slope (zero) over the full range.
    // ---------------------------------------------------------------------

    inline HygroThermFEM::SolidMaterialParams FiberglassBatts()
    {
        return {
            .name = "Fiberglass Batts",
            .thermalConductivityDry = 0.043,
            .density = 8.8,
            .porosity = 0.999,
            .heatCapacity = 840.0,
            .diffusionResistanceFactor = 1.21,
            .thermalConductivityMoistureDependent = {{0.0, 0.043}, {13.4, 0.043}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{10.0, 0.043}, {80.0, 0.043}},
            .temperatureDependentMeasurementHumidity = 0,
            // Both Suction and Redistribution curves are trivially empty in the XML
            .liquidTransportCurve = {{0, 0}},
            .sorptionCurve = {{0,     0},
                              {0.1,   0.002},
                              {0.5,   0.014},
                              {0.8,   0.53},
                              {0.9,   1.14},
                              {0.97,  3.37},
                              {0.99,  6.79},
                              {0.999, 12.2},
                              {1,     13.4}}
        };
    }

    inline HygroThermFEM::SolidMaterialParams GypsumBoardInterior()
    {
        return {
            .name = "Gypsum Board Interior",
            .thermalConductivityDry = 0.2,
            .density = 625.0,
            .porosity = 0.73,
            .heatCapacity = 850.0,
            .diffusionResistanceFactor = 8.33,
            .thermalConductivityMoistureDependent = {{0.0, 0.2}, {370.0, 0.2}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{10.0, 0.2}, {80.0, 0.2}},
            .temperatureDependentMeasurementHumidity = 0,
            // Suction curve from XML (Redistribution is identical)
            .liquidTransportCurve = {{0,   1.85e-10},
                                     {9.5, 1.85e-10},
                                     {370, 1.59e-07}},
            .sorptionCurve = {{0,    0},
                              {0.33, 5},
                              {0.75, 7},
                              {0.97, 18},
                              {1,    370}}
        };
    }

    inline HygroThermFEM::SolidMaterialParams LaminatedPanel()
    {
        return {
            .name = "Laminated panel",
            .thermalConductivityDry = 0.125,
            .density = 450.0,
            .porosity = 0.55,
            .heatCapacity = 1400.0,
            .diffusionResistanceFactor = 203.0,
            .thermalConductivityMoistureDependent = {{0.0, 0.12}, {534.0, 0.12}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{10.0, 0.12}, {80.0, 0.12}},
            .temperatureDependentMeasurementHumidity = 0,
            // Suction curve from XML
            .liquidTransportCurve = {{0,   0},
                                     {73,  4e-12},
                                     {534, 5e-12}},
            .sorptionCurve = {{0,      0},
                              {0.1,    37},
                              {0.3,    45},
                              {0.5,    53},
                              {0.7,    64},
                              {0.8,    73},
                              {0.9,    89},
                              {0.93,   98},
                              {0.95,   107},
                              {0.99,   159},
                              {0.995,  185},
                              {0.999,  251},
                              {0.9995, 281},
                              {0.9999, 347},
                              {1,      534}}
        };
    }

    inline HygroThermFEM::SolidMaterialParams Stucco()
    {
        return {
            .name = "Stucco",
            .thermalConductivityDry = 0.85,
            .density = 1800.0,
            .porosity = 0.24,
            .heatCapacity = 850.0,
            .diffusionResistanceFactor = 19.0,
            .thermalConductivityMoistureDependent = {{0.0, 0.8}, {210.0, 0.8}},
            .moistureDependentMeasurementTemperature = 0,
            .thermalConductivityTemperatureDependent = {{10.0, 0.8}, {80.0, 0.8}},
            .temperatureDependentMeasurementHumidity = 0,
            // Suction is empty in XML; fall back to Redistribution curve
            .liquidTransportCurve = {{0,   0},
                                     {20,  1e-10},
                                     {210, 7e-09}},
            .sorptionCurve = {{0,      0},
                              {0.5,    30},
                              {0.8,    45},
                              {0.9,    65},
                              {0.99,   95},
                              {0.999,  110},
                              {0.9995, 140},
                              {0.9999, 200},
                              {1,      210}}
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
