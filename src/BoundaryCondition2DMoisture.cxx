#include "BoundaryCondition2DMoisture.hxx"

namespace HygroThermFEM
{
    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////
    MoistureBCTARPHc::MoistureBCTARPHc(size_t index1,
                                       size_t index2,
                                       const std::string & materialName,
                                       const TARPCoefficients & varHCCoeff,
                                       double surfaceTilt) :
        IMoistureBC(index1,
                    index2,
                    materialName,
                    varHCCoeff.AirHumidity,
                    varHCCoeff.AirTemperature,
                    ConvectionModelFactory::createTARPFilmCoefficient(
                      m_Nodes, varHCCoeff.AirTemperature, surfaceTilt))
    {}

    /////////////////////////////////////////////////////
    /// MoistureBCASHRAEInside
    /////////////////////////////////////////////////////
    MoistureBCASHRAEInside::MoistureBCASHRAEInside(size_t index1,
                                                   size_t index2,
                                                   const std::string & materialName,
                                                   const ASHRAEInsideCoefficients & coeffs,
                                                   const double surfaceHeight,
                                                   const double surfaceTilt) :
        IMoistureBC(
          index1,
          index2,
          materialName,
          coeffs.AirHumidity,
          coeffs.AirTemperature,
          ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
            m_Nodes, coeffs.AirTemperature, surfaceTilt, surfaceHeight, coeffs.AirPressure))
    {}


    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////
    MoistureBCFixedHc::MoistureBCFixedHc(size_t index1,
                                         size_t index2,
                                         const std::string & materialName,
                                         const FixedBCHCCoefficients & fixedBchcCoefficients) :
        IMoistureBC(index1,
                    index2,
                    materialName,
                    fixedBchcCoefficients.AirHumidity,
                    fixedBchcCoefficients.AirTemperature,
                    ConvectionModelFactory::createFixedFilmCoefficient(
                      m_Nodes, fixedBchcCoefficients.ConvectionCoefficient))
    {}
}   // namespace HygroThermFEM
