#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"
#include "Materials.hxx"
#include "BoundaryConditionCoefficients.hxx"
#include "ConvectiveCoefficient.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// IConvectiveBCBase - Base class for convection-based BCs
    ////////////////////////////////////////////////////////

    //! \brief Base class for boundary conditions that use convective heat transfer models.
    //!
    //! Provides common members for air temperature, humidity, and convective coefficient
    //! calculator used by both thermal (IConvectionBC) and moisture (IMoistureBC) domains.
    class IConvectiveBCBase : public IBCLinear2D
    {
    public:
        virtual ~IConvectiveBCBase() = default;

    protected:
        //! Constructor for convective boundary condition base.
        IConvectiveBCBase(size_t index1,
                          size_t index2,
                          double airTemperature,
                          double airHumidity,
                          std::unique_ptr<IConvectiveCoefficient> convectiveCoeffCalc);

        const double m_AirTemperature;
        const double m_AirHumidity;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
    };

    ////////////////////////////////////////////////////////
    /// IConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Interface boundary condition class that is used to create various convection boundary
    //! condition models
    //!
    //! Simple convection boundary condition that is used in thermal module.
    class IConvectionBC : public IConvectiveBCBase
    {
    public:
        //! Constructor for convection boundary condition that is common between models
        IConvectionBC(size_t index1,
                      size_t index2,
                      double airTemperature,
                      std::unique_ptr<IConvectiveCoefficient> convectiveCoeffCalc,
                      double airHumidity = 0,
                      bool simulateMoisture = true);

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates H matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        //! Need to include energy from vapor flux energy if moisture is present.
        bool m_SimulateVaporFluxEnergy;
    };

    /////////////////////////////////////////////////////
    /// IMoistureBC
    /////////////////////////////////////////////////////

    //! \brief Base class for moisture boundary conditions. The only difference in inherited classes
    //! is how convection heat transfer coefficient is calculated.
    //!
    //! This class handles all matrices and vectors generation that is common for moisture bcs
    class IMoistureBC : public IConvectiveBCBase
    {
    public:
        //! Construction for moisture boundary condition.
        IMoistureBC(Materials & materialPool,
                    size_t index1,
                    size_t index2,
                    const std::string & materialName,
                    double airHumidity,
                    double airTemperature,
                    std::unique_ptr<IConvectiveCoefficient> convectiveCoeffCalc);

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        const IMaterial & m_Material;
    };

}   // namespace HygroThermFEM
