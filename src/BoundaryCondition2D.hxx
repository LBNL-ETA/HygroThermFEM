#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"
#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// IConvectiveCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Interface for convective coefficient calculations
    //!
    //! The most important functionality of this interface is that it keeps all the nodes
    //! associated with it.
    class IConvectiveCoefficient
    {
    public:
        explicit IConvectiveCoefficient(const INodes & nodes);
        virtual ~IConvectiveCoefficient() = default;
        IConvectiveCoefficient(const IConvectiveCoefficient & other) = default;
        IConvectiveCoefficient(IConvectiveCoefficient && other) = default;
        IConvectiveCoefficient & operator=(const IConvectiveCoefficient & other) = delete;
        IConvectiveCoefficient & operator=(IConvectiveCoefficient && other) = delete;
        [[nodiscard]] virtual std::vector<double> convectiveCoefficients() const = 0;

        //! \brief Water vapor transfer coefficient is always the same and it does not depend on
        //! how convective film coefficient is calculated
        //!
        //! \return Array that contains calculated values for each node associated with the boundary
        [[nodiscard]] std::vector<double> waterVaporTransferCoefficient() const;

    protected:
        //! Nodes associated with the boundary
        const INodes & m_Nodes;
    };

    ////////////////////////////////////////////////////////
    /// FixedConvectionCoefficient
    ////////////////////////////////////////////////////////

    class FixedConvectionCoefficient : public IConvectiveCoefficient
    {
    public:
        FixedConvectionCoefficient(const INodes & nodes, double convectionFilmCoefficient);

        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;

    private:
        double m_ConvectionFilmCoefficient;
    };

    ////////////////////////////////////////////////////////
    /// IVariableConvectiveCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Convective heat transfer coefficient is calculated based on ambient temperature
    class IVariableConvectiveCoefficient : public IConvectiveCoefficient
    {
    public:
        IVariableConvectiveCoefficient(const INodes & nodes,
                                       double airTemperature,
                                       double surfaceTilt = 0);

    protected:
        double m_AirTemperature{0};
        double m_SurfaceTilt{0};
    };

    ////////////////////////////////////////////////////////
    /// TARPFilmCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Film coefficient calculated based on comprehensive natural convection model (TARP)
    class TARPFilmCoefficient : public IVariableConvectiveCoefficient
    {
    public:
        //! \brief TARP algorithm requires air temperature and surface tilt.
        //!
        //! \param nodes Nodes associated with the boundary
        //! \param airTemperature Air temperature of the ambient [degrees Celsius]
        //! \param surfaceTilt Surface tilt for which film coefficient is being calculated [degrees]
        TARPFilmCoefficient(const INodes & nodes, double airTemperature, double surfaceTilt = 90);
        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;
    };

    ////////////////////////////////////////////////////////
    /// ASHRAEInsideFilmCoefficient
    ////////////////////////////////////////////////////////

    class ASHRAEInsideFilmCoefficient : public IVariableConvectiveCoefficient
    {
    public:
        ASHRAEInsideFilmCoefficient(const INodes & nodes,
                                    double airTemperature,
                                    double surfaceTilt,
                                    double mSurfaceHeight);

        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;

    private:
        double m_SurfaceHeight;
    };

    ////////////////////////////////////////////////////////
    /// ConvectionModelFactory
    ////////////////////////////////////////////////////////

    //! \brief Factory class for convection calculation model
    class ConvectionModelFactory
    {
    public:
        static std::unique_ptr<IConvectiveCoefficient>
          createConstantFilmCoefficient(const INodes & nodes, double filmCoefficient);

        static std::unique_ptr<IConvectiveCoefficient> createTARPFilmCoefficient(
          const INodes & nodes, double airTemperature, double surfaceTilt = 90);
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

        //! If moisture is simulated then calculate energy from moisture flux. Otherwise, don't do
        //! it.
        bool m_SimulateMoisture;
    };

    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    class ConstantConvectionBC : public IConvectionBC
    {
    public:
        ConstantConvectionBC(size_t index1,
                             size_t index2,
                             const FixedBCHCCoefficients & fixedBCHCCoefficients,
                             bool t_CalculateMoisture = true);
    };

    ////////////////////////////////////////////////////////
    /// TARPConvectionBC
    ////////////////////////////////////////////////////////
    class TARPConvectionBC : public IConvectionBC
    {
    public:
        TARPConvectionBC(size_t index1,
                         size_t index2,
                         const VariableBCHCCoefficients & varHCCoeff,
                         bool t_CalculateMoisture = true);
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
        SquareMatrix H_Matrix() const override;

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
          double t_RadiationTemperature   //!< Outside radiation temperature in celsius.
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        SquareMatrix H_Matrix() const override;

    private:
        //! Radiative convective coefficient that needs to be calculated based on current
        //! temperatures
        std::vector<double> HRadiative() const;

        double m_RadiationTemperature;
        double m_Emissivity;
    };

    ///////////////////////////////////////////////////////
    /// SimplifiedRadiationBC
    ///////////////////////////////////////////////////////
    class LinearizedRadiationBC : public IBCLinear2D
    {
    public:
        LinearizedRadiationBC(size_t index1,
                              size_t index2,
                              const LinearizedRadiationBCCoefficients & linearRadBC);

        std::vector<double> R_Vector() const override;

        SquareMatrix H_Matrix() const override;

    private:
        double m_RadiationCoefficient;
        double m_RadiationTemperature;
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
        SquareMatrix H_Matrix() const override;

    protected:
        double m_AirHumidity;
        double m_AirTemperature;
        const IMaterial & m_Material;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
    };

    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient.
    class MoistureBCTARPHc : public IMoistureBC
    {
    public:
        MoistureBCTARPHc(size_t index1,
                         size_t index2,
                         const std::string & materialName,
                         const VariableBCHCCoefficients & varHCCoeff);
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
                          const FixedBCHCCoefficients & fixedBchcCoefficients);
    };

}   // namespace HygroThermFEM
