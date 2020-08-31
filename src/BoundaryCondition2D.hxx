#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"
#include "BoundaryConditionCoefficients.hxx"
#include "ConvectiveCoefficient.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// IConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Interface boundary condition class that is used to create various convection boundary
    //! condition models
    //!
    //! Simple convection boundary condition that is used in thermal module.
    class IConvectionBC : public IBCLinear2D
    {
    public:
        //! Constructor for convection boundary condition that is common between models
        IConvectionBC(size_t index1,
                      size_t index2,
                      double t_AirTemperature,
                      std::unique_ptr<IConvectiveCoefficient>,
                      double t_AirHumidity = 0,
                      bool t_SimulateMoisture = true);

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates H matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        const double m_AirTemperature;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
        const double m_AirHumidity;

        //! Need to include energy from vapor flux energy if moisture is present.
        bool m_SimulateVaporFluxEnergy;
    };    

    /////////////////////////////////////////////////////
    /// MoistureBC
    /////////////////////////////////////////////////////

    //! \brief Base class for moisture boundary conditions. The only difference in inherited classes
    //! is how convection heat transfer coefficient is calculated.
    //!
    //! This class handles all matrices and vectors generation that is common for moisture bcs
    class IMoistureBC : public IBCLinear2D
    {
    public:
        //! Construction for moisture boundary condition.
        IMoistureBC(size_t index1,
                    size_t index2,
                    const std::string & materialName,
                    double t_AirHumidity,
                    double t_AirTemperature,
                    std::unique_ptr<IConvectiveCoefficient> model);

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        double m_AirHumidity;
        double m_AirTemperature;
        const IMaterial & m_Material;

        //! \brief Convective coefficients calculations can be performed with different algorithms
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
    };    

}   // namespace HygroThermFEM
