#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"

namespace HygroThermFEM
{
    enum class ConvectionModel
    {
        Fixed,
        Variable
    };

    ////////////////////////////////////////////////////////
    /// IConvectiveCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Interface for convective coefficient calculations
    class IConvectiveCoefficient
    {
    public:
        IConvectiveCoefficient(const INodes & nodes, double ambientVariable);
        virtual ~IConvectiveCoefficient() = default;
        IConvectiveCoefficient(const IConvectiveCoefficient & other) = default;
        IConvectiveCoefficient(IConvectiveCoefficient && other) = default;
        IConvectiveCoefficient & operator=(const IConvectiveCoefficient & other) = delete;
        IConvectiveCoefficient & operator=(IConvectiveCoefficient && other) = delete;
        virtual std::vector<double> convectiveCoefficients() const = 0;
        std::vector<double> betaConv() const;

    protected:
        const INodes & m_Nodes;
        const double m_AmbientVariable;
    };

    ////////////////////////////////////////////////////////
    /// VariableConvectionCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Convective heat transfer coefficient is calculated based on ambient temperature
    class VariableConvectionCoefficient : public IConvectiveCoefficient
    {
    public:
        VariableConvectionCoefficient(const INodes & nodes, double ambientVariable);

        std::vector<double> convectiveCoefficients() const override;
    };

    ////////////////////////////////////////////////////////
    /// FixedConvectionCoefficient
    ////////////////////////////////////////////////////////

    class FixedConvectionCoefficient : public IConvectiveCoefficient
    {
    public:
        FixedConvectionCoefficient(const INodes & nodes, double ambientVariable);

        std::vector<double> convectiveCoefficients() const override;
    };

    ////////////////////////////////////////////////////////
    /// ConvectionModelFactory
    ////////////////////////////////////////////////////////

    //! \brief Factory class for convection calculation model
    class ConvectionModelFactory
    {
    public:
        static std::unique_ptr<IConvectiveCoefficient>
          create(ConvectionModel model, const INodes & nodes, double ambientVariable);
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
        IConvectionBC(size_t index1,
                      size_t index2,
                      double t_AirTemperature,
                      std::unique_ptr<IConvectiveCoefficient>,
                      double t_AirHumidity = 0);

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    protected:
        const double m_AirTemperature;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
        const double m_AirHumidity;
    };

    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    class ConstantConvectionBC : public IConvectionBC
    {
    public:
        ConstantConvectionBC(size_t index1,
                             size_t index2,
                             double t_AirTemperature,
                             double t_ConvectionCoefficient,
                             double t_AirHumidity = 0);
    };

    ////////////////////////////////////////////////////////
    /// VariableConvectionBC
    ////////////////////////////////////////////////////////
    class VariableConvectionBC : public IConvectionBC
    {
    public:
        VariableConvectionBC(size_t index1,
                             size_t index2,
                             double t_AirTemperature,
                             double t_AirHumidity);
    };

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    //! Constant temperature boundary condition
    //!
    //! TemperatureBC will be just special case of convection BC with huge convectiveCoefficients
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
                    std::unique_ptr<IConvectiveCoefficient> model);

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    protected:

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

    private:
        const double m_ConvectiveCoefficient;
    };

}   // namespace HygroThermFEM
