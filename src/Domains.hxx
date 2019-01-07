#pragma once

#include <memory>
#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "Elements2D.hxx"

namespace MoisThermFEM
{
    //! \brief Interface that will keep all elements and boundary conditions together.
    //!
    //! One domain will solve single differential equation and therefore, single domain will
    //! represent thermal, moisture or pressure separately.
    class IDomain
    {
    public:
        //! Domain construction. It is necessary to set up base variable that will be considered
        //! unknown.
        explicit IDomain(
          const BaseVariable property   //!< State variable which will be considered unknown.
        );

        //! Calculates steady state for given data
        std::vector<double> steadyState();

        //! Calculates next timestep values from current (initial) values
        std::vector<double> transient(
          const std::vector<double> &
            currentStateValues,   //!< Current values of state variable or initial condition
          double t_DTime          //!< Timestep in transient solution
        );

        //! Returns flux in x and y direction
        std::vector< NodeFlux > flux() const;

        //! Adds element into domain
        virtual void createElement(
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< Material that will be assigned to the element
          ) = 0;

    protected:
        friend class MultiDomain;

        //! Helper class that will update property at nodes from new timestep values.
        void updateNodeValues(const std::vector<double> & values, const BaseVariable property);

        //! Forms left hand side matrix in steady state solution.
        FenestrationCommon::SquareMatrix steadyStateLeftHandSide();

        //! Form right hand side vector in stead state solution.
        std::vector<double> steadyStateRightHandSide() const;

        //! Forms mass, conductance and H (from boundary condition) matrices.
        FenestrationCommon::SquareMatrix transientM_K_H_Matrix(const double t_DTime);

        //! This function retrieves M*U+R vector (where U is state variable)
        std::vector<double> transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                                 const double t_DTime);

        //! Returns if domain problem is linear.
        bool isLinear() const;

        //! Some domains require post-processing of results. Good example is
        //! moisture domain where humidity cannot go over 1.0 or lower than one.
        //! With certain set of boundary conditions and long enough time-step,
        //! solution can achieve such state and post processing should prevent it.
        virtual void postProcess(std::vector<double> & solution) const;

        BaseVariable m_Property;
        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;
    };

    //! \brief Domain class for solving temperature solution.
    class ThermalDomain : public IDomain
    {
    public:
        //! Simple constructor
        ThermalDomain();

        //! Creation of convection boundary condition
        void createConvectionBC(
          size_t index1,                    //!< Node 1 index
          size_t index2,                    //!< Node 2 index
          double t_ConvectionCoefficient,   //!< Heat transfer convection coefficient
          double t_AirTemperature           //!< Outside air temperature
        );

        //! Creation of temperature boundary condition
        void createTemperatureBC(size_t index1,    //!< Node 1 index
                                 size_t index2,    //!< Node 2 index
                                 double t_Temp1,   //!< Constant temperature at node 1
                                 double t_Temp2    //!< Constant temperature at node 2
        );

        //! Creation of temperature boundary condition
        void createTemperatureBC(size_t index1,   //!< Node 1 index
                                 size_t index2,   //!< Node 2 index
                                 double t_Temp    //!< Constant temperature in both nodes.
        );

        //! Creation of flux boundary condition
        void createFluxBC(size_t index1,   //!< Node 1 index
                          size_t index2,   //!< Node 2 index
                          double t_Flux    //!< Constant flux in both nodes.
        );

        //! Creation of black body radiation
        void createBlackBodyRadiationBC(size_t index1,         //!< Node 1 index
                                        size_t index2,         //!< Node 2 index
                                        double t_Emissivity,   //! Emissivity at boundary condition
                                        double t_RadiationTemperature   //! Radiation temperature
        );

        //! Creates and adds element into domain.
        virtual void createElement(size_t index1,                     //!< Node 1 index
                                   size_t index2,                     //!< Node 2 index
                                   size_t index3,                     //!< Node 3 index
                                   size_t index4,                     //!< Node 4 index
                                   const std::string & materialName   //!< Material name
                                   ) override;
    };

    //! \brief Domain class for solving humidity distribution.
    class MoistureDomain : public IDomain
    {
    public:
        //! Simple constructor
        MoistureDomain();

        //! Creates moisture boundary conditions.
		void createMoistureBC( const size_t index1, const size_t index2,
							   const double t_AirHumidity, const double t_AirTemperature );

        //! Creates and adds element into domain.
        virtual void createElement(size_t index1,                     //!< Node 1 index
                                   size_t index2,                     //!< Node 2 index
                                   size_t index3,                     //!< Node 3 index
                                   size_t index4,                     //!< Node 4 index
                                   const std::string & materialName   //!< Material name
                                   ) override;
	protected:
		void postProcess( std::vector< double > & solution ) const override;
	};

}   // namespace MoisThermFEM