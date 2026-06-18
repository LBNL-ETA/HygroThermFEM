#pragma once

#include <memory>
#include <vector>
#include <map>

namespace HygroThermFEM
{
    class INodes;
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
    /// TARPFilmCoefficient
    ////////////////////////////////////////////////////////

    //! \brief Film coefficient calculated based on comprehensive natural convection model (TARP)
    class TARPFilmCoefficient : public IConvectiveCoefficient
    {
    public:
        //! \brief TARP algorithm requires air temperature and surface tilt.
        //!
        //! \param nodes Nodes associated with the boundary
        //! \param airTemperature Air temperature of the ambient [degrees Celsius]
        //! \param surfaceTilt Surface tilt for which film coefficient is being calculated [degrees]
        TARPFilmCoefficient(const INodes & nodes, double airTemperature, double surfaceTilt = 90);
        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;

    private:
        double m_AirTemperature;
        double m_SurfaceTilt;
    };

    ////////////////////////////////////////////////////////
    /// ASHRAEInsideFilmCoefficient
    ////////////////////////////////////////////////////////

    class ASHRAEInsideFilmCoefficient : public IConvectiveCoefficient
    {
    public:
        ASHRAEInsideFilmCoefficient(const INodes & nodes,
                                    double airTemperature,
                                    double airPressure,
                                    double surfaceTilt,
                                    double surfaceHeight);

        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;

    private:
        double m_AirTemperature;
        double m_AirPressure;
        double m_SurfaceTilt;
        double m_SurfaceHeight;
    };

    ////////////////////////////////////////////////////////
    /// ASHRAEOutsideFilmCoefficient
    ////////////////////////////////////////////////////////

    class ASHRAEOutsideFilmCoefficient : public IConvectiveCoefficient
    {
    public:
        ASHRAEOutsideFilmCoefficient(const INodes & nodes, double windSpeed);

        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;

    private:
        double m_WindSpeed;
    };

    enum class WindDirection
    {
        Windward,
        Leeward
    };

    ////////////////////////////////////////////////////////
    /// YazdanianKlemsFilmCoefficient
    ////////////////////////////////////////////////////////

    class YazdanianKlemsFilmCoefficient : public IConvectiveCoefficient
    {
    public:
        YazdanianKlemsFilmCoefficient(const INodes & nodes,
                                      double airTemperature,
                                      double windSpeed,
                                      WindDirection direction);

        [[nodiscard]] std::vector<double> convectiveCoefficients() const override;

    private:
        double m_AirTemperature;
        double m_WindSpeed;
        WindDirection m_Direction;

        struct Coeffs
        {
            double A;
            double B;
        };

        std::map<WindDirection, Coeffs> coeffs{{WindDirection::Windward, {2.38, 0.89}},
                                               {WindDirection::Leeward, {2.86, 0.617}}};
    };

    ////////////////////////////////////////////////////////
    /// KimuraFilmCoefficient
    ////////////////////////////////////////////////////////

    class KimuraFilmCoefficient : public IConvectiveCoefficient
    {
    public:
        KimuraFilmCoefficient(const INodes & nodes, double mWindSpeed, WindDirection mDirection);

        std::vector<double> convectiveCoefficients() const override;

    private:
        double m_WindSpeed;
        WindDirection m_Direction;
    };

    ////////////////////////////////////////////////////////
    /// ConvectionModelFactory
    ////////////////////////////////////////////////////////

    //! \brief Factory class for convection calculation model
    class ConvectionModelFactory
    {
    public:
        static std::unique_ptr<IConvectiveCoefficient>
          createFixedFilmCoefficient(const INodes & nodes, double filmCoefficient);

        static std::unique_ptr<IConvectiveCoefficient> createTARPFilmCoefficient(
          const INodes & nodes, double airTemperature, double surfaceTilt = 90);

        static std::unique_ptr<IConvectiveCoefficient>
          createASHRAEOutsideFilmCoefficient(const INodes & nodes, double windSpeed);

        static std::unique_ptr<IConvectiveCoefficient>
          createASHRAEInsideFilmCoefficient(const INodes & nodes,
                                            double airTemperature,
                                            double surfaceTilt,
                                            double surfaceHeight,
                                            double airPressure = 101325);

        static std::unique_ptr<IConvectiveCoefficient> createYazdanianKlemsFilmCoefficient(
          const INodes & nodes, double airTemperature, double windSpeed, WindDirection direction);

        static std::unique_ptr<IConvectiveCoefficient> createKimuraFilmCoefficient(
          const INodes & nodes, double windSpeed, WindDirection direction);
    };
}