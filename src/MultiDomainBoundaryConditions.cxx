#include "MultiDomainBoundaryConditions.hxx"
#include "MultiDomain.hxx"

namespace HygroThermFEM
{
    void createBC_FixedHc(MultiDomain & domain,
                          const size_t index1,
                          const size_t index2,
                          const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        domain.thermalDomain.createBC_FixedHc(
          index1, index2, fixedBchcCoefficients, domain.simulateMoisture);
        domain.moistureDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients);
    }

    void createBC_FixedHc(MultiDomain & domain,
                          size_t index1,
                          size_t index2,
                          const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        domain.thermalDomain.createBC_FixedHc(
          index1, index2, fixedBchcCoefficients, domain.simulateMoisture);
        domain.moistureDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients);
    }

    void createBC_TARPHc(MultiDomain & domain,
                         size_t index1,
                         size_t index2,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt)
    {
        domain.thermalDomain.createBC_TARPHc(
          index1, index2, varHCCoeff, surfaceTilt, domain.simulateMoisture);

        domain.moistureDomain.createBC_TARPHc(index1, index2, varHCCoeff, surfaceTilt);
    }

    void createBC_TARPHc(MultiDomain & domain,
                         size_t index1,
                         size_t index2,
                         const std::vector<TARPCoefficients> & varHCCoeff,
                         double surfaceTilt)
    {
        domain.thermalDomain.createBC_TARPHc(
          index1, index2, varHCCoeff, surfaceTilt, domain.simulateMoisture);

        domain.moistureDomain.createBC_TARPHc(index1, index2, varHCCoeff, surfaceTilt);
    }

    void createBC_ASHRAEInsideHc(MultiDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt)
    {
        domain.thermalDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt, domain.simulateMoisture);
        domain.moistureDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt);
    }

    void createBC_ASHRAEInsideHc(MultiDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const std::vector<ASHRAEInsideCoefficients> & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt)
    {
        domain.thermalDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt, domain.simulateMoisture);
        domain.moistureDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt);
    }

    void createBC_ASHRAEOutsideHc(MultiDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff)
    {
        domain.thermalDomain.createBC_ASHRAEOutsideHc(
          index1, index2, coeff, domain.simulateMoisture);
        domain.moistureDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff);
    }

    void createBC_ASHRAEOutsideHc(MultiDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<ASHRAEOutsideCoefficients> & coeff)
    {
        domain.thermalDomain.createBC_ASHRAEOutsideHc(
          index1, index2, coeff, domain.simulateMoisture);
        domain.moistureDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff);
    }

    void createBC_YazdanianKlemsHc(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff)
    {
        domain.thermalDomain.createBC_YazdanianKlemsHc(
          index1, index2, coeff, domain.simulateMoisture);
        domain.moistureDomain.createBC_YazdanianKlemsHc(index1, index2, coeff);
    }

    void createBC_YazdanianKlemsHc(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<YazdanianKlemsCoefficients> & coeff)
    {
        domain.thermalDomain.createBC_YazdanianKlemsHc(
          index1, index2, coeff, domain.simulateMoisture);
        domain.moistureDomain.createBC_YazdanianKlemsHc(index1, index2, coeff);
    }

    void createBC_KimuraHc(MultiDomain & domain,
                           size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff)
    {
        domain.thermalDomain.createBC_KimuraHc(index1, index2, coeff, domain.simulateMoisture);
        domain.moistureDomain.createBC_KimuraHc(index1, index2, coeff);
    }

    void createBC_KimuraHc(MultiDomain & domain,
                           size_t index1,
                           size_t index2,
                           const std::vector<KimuraCoefficients> & coeff)
    {
        domain.thermalDomain.createBC_KimuraHc(index1, index2, coeff, domain.simulateMoisture);
        domain.moistureDomain.createBC_KimuraHc(index1, index2, coeff);
    }

    void createBC_FixedTemperature(MultiDomain & domain,
                                   const size_t index1,
                                   const size_t index2,
                                   const double t_Temp1,
                                   const double t_Temp2)
    {
        domain.thermalDomain.createBC_FixedTemperature(index1, index2, t_Temp1, t_Temp2);
    }

    void createBC_FixedTemperature(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<ConstantBCTemperatures> & temp)
    {
        domain.thermalDomain.createBC_FixedTemperature(index1, index2, temp);
    }

    void createBC_FixedTemperature(MultiDomain & domain,
                                   const size_t index1,
                                   const size_t index2,
                                   const double t_Temp)
    {
        domain.thermalDomain.createBC_FixedTemperature(index1, index2, t_Temp);
    }

    void createBC_FixedTemperature(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   std::vector<double> temp)
    {
        domain.thermalDomain.createBC_FixedTemperature(index1, index2, std::move(temp));
    }

    void createBC_FixedTemperatureAndHumidity(MultiDomain & domain,
                                              size_t index1,
                                              size_t index2,
                                              const TemperatureAndHumidity & values)
    {
        domain.thermalDomain.createBC_FixedTemperature(index1, index2, values.Temperature);
        domain.moistureDomain.createBC_FixedHumidity(index1, index2, values);
    }

    void createBC_FixedTemperatureAndHumidity(MultiDomain & domain,
                                              size_t index1,
                                              size_t index2,
                                              const std::vector<TemperatureAndHumidity> & values)
    {
        std::vector<double> temperatures(values.size());
        for(size_t i = 0u; i < values.size(); ++i)
        {
            temperatures[i] = values[i].Temperature;
        }
        domain.thermalDomain.createBC_FixedTemperature(index1, index2, temperatures);
        domain.moistureDomain.createBC_FixedHumidity(index1, index2, values);
    }

    void createBC_FixedHeatFlux(MultiDomain & domain, size_t index1, size_t index2, double t_Flux)
    {
        domain.thermalDomain.createBC_FixedFlux(index1, index2, t_Flux);
    }

    void createBC_FixedHeatFlux(MultiDomain & domain,
                                size_t index1,
                                size_t index2,
                                const std::vector<double> & t_Flux)
    {
        domain.thermalDomain.createBC_FixedFlux(index1, index2, t_Flux);
    }

    void createBC_BlackBodyRadiation(MultiDomain & domain,
                                     size_t index1,
                                     size_t index2,
                                     double t_Emissivity,
                                     double t_RadiationTemperature)
    {
        domain.thermalDomain.createBC_BlackBodyRadiation(
          index1, index2, t_Emissivity, t_RadiationTemperature);
    }

    void
      createBC_BlackBodyRadiation(MultiDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        domain.thermalDomain.createBC_BlackBodyRadiation(index1, index2, radCoeffs);
    }

    void createBC_LinearizedRadiation(MultiDomain & domain,
                                      const size_t index1,
                                      const size_t index2,
                                      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        domain.thermalDomain.createBC_LinearizedRadiation(index1, index2, linearRadBC);
    }

    void createBC_LinearizedRadiation(
      MultiDomain & domain,
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC)
    {
        domain.thermalDomain.createBC_LinearizedRadiation(index1, index2, linearRadBC);
    }
}   // namespace HygroThermFEM