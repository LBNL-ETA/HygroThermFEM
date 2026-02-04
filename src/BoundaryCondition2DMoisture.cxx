#include "BoundaryCondition2DMoisture.hxx"

namespace HygroThermFEM
{
    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////
    MoistureBCTARPHc::MoistureBCTARPHc(Nodes & nodePool,
                                       Materials & materialPool,
                                       const size_t index1,
                                       const size_t index2,
                                       const std::string & materialName,
                                       const TARPCoefficients & varHCCoeff,
                                       const double surfaceTilt) :
        IMoistureBC(nodePool,
                    materialPool,
                    index1,
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
    MoistureBCASHRAEInside::MoistureBCASHRAEInside(Nodes & nodePool,
                                                   Materials & materialPool,
                                                   const size_t index1,
                                                   const size_t index2,
                                                   const std::string & materialName,
                                                   const ASHRAEInsideCoefficients & coeffs,
                                                   const double surfaceHeight,
                                                   const double surfaceTilt) :
        IMoistureBC(
          nodePool,
          materialPool,
          index1,
          index2,
          materialName,
          coeffs.AirHumidity,
          coeffs.AirTemperature,
          ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
            m_Nodes, coeffs.AirTemperature, surfaceTilt, surfaceHeight, coeffs.AirPressure))
    {}

    /////////////////////////////////////////////////////
    /// MoistureBCASHRAEOutside
    /////////////////////////////////////////////////////
    MoistureBCASHRAEOutside::MoistureBCASHRAEOutside(Nodes & nodePool,
                                                     Materials & materialPool,
                                                     const size_t index1,
                                                     const size_t index2,
                                                     const std::string & materialName,
                                                     const ASHRAEOutsideCoefficients & coeffs) :
        IMoistureBC(
          nodePool,
          materialPool,
          index1,
          index2,
          materialName,
          coeffs.AirHumidity,
          coeffs.AirTemperature,
          ConvectionModelFactory::createASHRAEOutsideFilmCoefficient(m_Nodes, coeffs.WindSpeed))
    {}

    /////////////////////////////////////////////////////
    /// MoistureYazdanianKlemsBC
    /////////////////////////////////////////////////////
    MoistureYazdanianKlemsBC::MoistureYazdanianKlemsBC(Nodes & nodePool,
                                                       Materials & materialPool,
                                                       const size_t index1,
                                                       const size_t index2,
                                                       const std::string & materialName,
                                                       const YazdanianKlemsCoefficients & coeffs) :
        IMoistureBC(nodePool,
                    materialPool,
                    index1,
                    index2,
                    materialName,
                    coeffs.AirHumidity,
                    coeffs.AirTemperature,
                    ConvectionModelFactory::createYazdanianKlemsFilmCoefficient(
                      m_Nodes, coeffs.AirTemperature, coeffs.WindSpeed, coeffs.WindDir))
    {}

    /////////////////////////////////////////////////////
    /// MoistureKimuraBC
    /////////////////////////////////////////////////////
    MoistureKimuraBC::MoistureKimuraBC(Nodes & nodePool,
                                       Materials & materialPool,
                                       const size_t index1,
                                       const size_t index2,
                                       const std::string & materialName,
                                       const KimuraCoefficients & coeffs) :
        IMoistureBC(nodePool,
                    materialPool,
                    index1,
                    index2,
                    materialName,
                    coeffs.AirHumidity,
                    coeffs.AirTemperature,
                    ConvectionModelFactory::createKimuraFilmCoefficient(
                      m_Nodes, coeffs.WindSpeed, coeffs.WindDir))
    {}

    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////
    MoistureBCFixedHc::MoistureBCFixedHc(Nodes & nodePool,
                                         Materials & materialPool,
                                         const size_t index1,
                                         const size_t index2,
                                         const std::string & materialName,
                                         const FixedBCHCCoefficients & fixedBchcCoefficients) :
        IMoistureBC(nodePool,
                    materialPool,
                    index1,
                    index2,
                    materialName,
                    fixedBchcCoefficients.AirHumidity,
                    fixedBchcCoefficients.AirTemperature,
                    ConvectionModelFactory::createFixedFilmCoefficient(
                      m_Nodes, fixedBchcCoefficients.ConvectionCoefficient))
    {}
}   // namespace HygroThermFEM
