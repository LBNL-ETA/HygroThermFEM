#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"

namespace MoisThermFEM
{
    enum class ConvectionModel {Fixed, Variable};

    //! \brief Interface for convective coefficient calculations
    class IConvectiveCoefficient
    {
    public:
        virtual ~IConvectiveCoefficient() = default;
        virtual std::vector<double> value(const INodes & nodes, const double variable) const = 0;
    };

    ////////////////////////////////////////////////////////
    /// VariableConvectionCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Convective heat transfer coefficient is calculated based on ambient temperature
    class VariableConvectionCoefficient : public IConvectiveCoefficient
    {
    public:
        VariableConvectionCoefficient();
        std::vector<double> value(const INodes & nodes, const double ambientTemperature) const override;
    };

    ////////////////////////////////////////////////////////
    /// FixedConvectionCoefficient
    ////////////////////////////////////////////////////////

    class FixedConvectionCoefficient : public IConvectiveCoefficient
    {
    public:
        FixedConvectionCoefficient();
        std::vector<double> value(const INodes & nodes, const double convectiveCoefficient) const override;
    };

    ////////////////////////////////////////////////////////
    /// ConvectionModelFactory
    ////////////////////////////////////////////////////////

    //! \brief Factory class for convection calculation model
    class ConvectionModelFactory
    {
    public:
        static std::unique_ptr<IConvectiveCoefficient> create(const ConvectionModel model);
    };

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
        IConvectionBC(size_t index1, size_t index2, double t_AirTemperature, ConvectionModel model);

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    protected:
        const double m_AirTemperature;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
        virtual std::vector<double> convectionCoefficients() const = 0;
    };

    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    class ConstantConvectionBC : public IConvectionBC
    {
    public:
        ConstantConvectionBC(size_t index1, size_t index2,
                                     double t_AirTemperature,
                                     const double m_ConvectionCoefficient);
    protected:
        std::vector<double> convectionCoefficients() const override;
    private:
        const double m_ConvectionCoefficient;
    };

    ////////////////////////////////////////////////////////
    /// VariableConvectionBC
    ////////////////////////////////////////////////////////
    class VariableConvectionBC : public IConvectionBC
    {
    public:
        VariableConvectionBC(size_t index1, size_t index2, double t_AirTemperature);
    protected:
        std::vector<double> convectionCoefficients() const override;
    };

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    //! Constant temperature boundary condition
    //!
    //! TemperatureBC will be just special case of convection BC with huge value
    //! for film coefficients. It is used only in thermal module.
    class TemperatureBC : public ConstantConvectionBC
    {
    public:
        //! Construction of temperature boundary condition with two nodes and identical temperatures
        //! in both of them
        TemperatureBC(size_t index1,              //!< Node 1 index
                      size_t index2,              //!< Node 2 index
                      double t_NodeTemperatures   //!< Temperature at both nodes
        );

        //! Construction of temperature boundary condition with two nodes and temperature at each of
        //! them
        TemperatureBC(size_t index1,    //!< Node 1 index
                      size_t index2,    //!< Node 2 index
                      double t_Temp1,   //!< Temperature at node 1
                      double t_Temp2    //!< Temperature at node 2
        );
    };

    ////////////////////////////////////////////////////////
    /// Flux BC
    ////////////////////////////////////////////////////////

    //! \brief Constant flux boundary condition.
    //!
    //! Boundary condition is used only in thermal module.
    class FluxBC : public IBCLinear2D
    {
    public:
        //! Construction of constant flux boundary condition with two nodes and constant flux.
        FluxBC(size_t index1,   //!< Node 1 index
               size_t index2,   //!< Node 2 index
               double t_Flux    //!< Constant flux at both nodes
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    private:
        double m_Flux;
    };

    ///////////////////////////////////////////////////////
    /// BlackBodyRadiationBC
    ///////////////////////////////////////////////////////

    //! \brief Simple black body radiation boundary condition.
    //!
    //! Boundary condition is used only in thermal module.
    class BlackBodyRadiationBC : public IBCLinear2D
    {
    public:
        //! Black body radiation boundary condition with two nodes, emissivity and radiation
        //! temperature
        BlackBodyRadiationBC(
          size_t index1,                  //!< Node 1 index
          size_t index2,                  //!< Node 2 index
          double t_Emissivity,            //!< Boundary condition surface emissivity at both nodes.
          double t_RadiationTemperature   //!< Outside radiation temperature.
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    private:
        //! Radiative convective coefficient that needs to be calculated based on current
        //! temperatures
        std::vector<double> HRadiative() const;

        double m_RadiationTemperature;
        double m_Emissivity;
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
                    ConvectionModel model);

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    protected:
        //! Calculates beta convective exterior (see technical document)
        std::vector<double> betaConv() const;
        virtual std::vector<double> convectiveCoefficient() const = 0;

        double m_AirHumidity;
        double m_AirTemperature;
        const Material & m_Material;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
    };

    /////////////////////////////////////////////////////
    /// MoistureBCCalculatedHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient.
    class MoistureBCVariableHc : public IMoistureBC
    {
    public:
        MoistureBCVariableHc(size_t index1,
                               size_t index2,
                               const std::string & materialName,
                               double t_AirHumidity,
                               double t_AirTemperature);
    protected:
        std::vector<double> convectiveCoefficient() const override;
    };

    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition with fixed convective coefficient
    class MoistureBCFixedHc : public IMoistureBC
    {
    public:
        MoistureBCFixedHc(size_t index1,
                          size_t index2,
                          const std::string & materialName,
                          double t_AirHumidity,
                          double t_AirTemperature,
                          double m_ConvectiveCoefficient);
    protected:
        std::vector<double> convectiveCoefficient() const override;

    private:
        const double m_ConvectiveCoefficient;
    };

}   // namespace MoisThermFEM
