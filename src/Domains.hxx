#pragma once

#include <memory>
#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "Elements2D.hxx"

namespace MoisThermFEM
{
    /// Class to keep elements and boundary conditions together. One domain per Thermal, Mass and
    /// Pressure governing equations.
    class Domain
    {
    public:
        explicit Domain(const BaseVariable property);

        /// Calculates steady state solution
        std::vector<double> steadyState();

        /// Calculates next timestep value from current values
        std::vector<double> transient(const std::vector<double> & currentStateValues,
                                      double t_DTime);

        virtual void createElement(const size_t index1,
                                   const size_t index2,
                                   const size_t index3,
                                   const size_t index4,
                                   const std::string & materialName) = 0;

        IElementLinear2D * findElement( const size_t index1, const size_t index2 );

    protected:
        friend class MultiDomain;

        void updateNodeValues(const std::vector<double> & values, const BaseVariable property);

        FenestrationCommon::SquareMatrix steadyStateLeftHandSide();
        std::vector<double> steadyStateRightHandSide() const;

        /// In matrix equations some structures are showing up in both (linear and nonlinear) cases
        /// and those matrix operations are separated into functions.
        /// This function retrieves M+K+H matrix
        FenestrationCommon::SquareMatrix transientM_K_H_Matrix(const double t_DTime);
        /// FenestrationCommon::SquareMatrix< double > transientDH_Matrix();

        /// This function retrieves M*U+R vector (where U is state variable)
        std::vector<double> transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                                 const double t_DTime);

        bool isLinear() const;

        BaseVariable m_Property;
        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;
    };

    class ThermalDomain : public Domain
    {
    public:
        ThermalDomain();

        void createConvectionBC(const size_t index1,
                                const size_t index2,
                                double t_ConvectionCoefficient,
                                double t_AirTemperature);

        void createTemperatureBC(const size_t index1,
                                 const size_t index2,
                                 double t_Temp1,
                                 double t_Temp2);

        void createTemperatureBC(const size_t index1, const size_t index2, const double t_Temp);

        void createFluxBC(const size_t index1, const size_t index2, const double t_Flux);

        void createBlackBodyRadiationBC(const size_t index1,
                                        const size_t index2,
                                        const double t_Emissivity,
                                        const double t_RadiationTemperature);

        virtual void createElement(const size_t index1,
                                   const size_t index2,
                                   const size_t index3,
                                   const size_t index4,
                                   const std::string & materialName) override;
    };

    class MoistureDomain : public Domain
    {
    public:
        MoistureDomain();

        void createMoistureBC(const size_t index1,
                              const size_t index2,
                              const double t_ConvectiveCoefficient,
                              const double t_AirHumidity,
                              const double t_AirTemperature);

        virtual void createElement(const size_t index1,
                                   const size_t index2,
                                   const size_t index3,
                                   const size_t index4,
                                   const std::string & materialName) override;
    };

}   // namespace MoisThermFEM