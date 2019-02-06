#pragma once

#include <vector>

namespace HygroThermFEM
{
    struct LocalPoint2D;
    struct LocalPoint1D;

    //////////////////////////////////////////////////////////////////////////////////////////////
    ////  ILocalShapeFunctions2DQuadrilateral
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Interface for quadrilater shape function.
    //!
    //! Interface for quadrilateral shape functions and its derivatives for point in two
    //! dimensional local coordinate system range from (-1, 1) in both directions
    class ILocalShapeFunctions2DQuadrilateral
    {
    public:
        virtual ~ILocalShapeFunctions2DQuadrilateral() = default;

        explicit ILocalShapeFunctions2DQuadrilateral(const LocalPoint2D & t_Point);

        virtual double Psi(size_t Index) const final;

        virtual double PsiDKsi(size_t Index) const final;
        virtual double PsiDEta(size_t Index) const final;

        virtual const std::vector< double > & VPsi() const final;
        virtual const std::vector< double > & VPsiDKsi() const final;
        virtual const std::vector< double > & VPsiDEta() const final;

    protected:
        double m_Ksi;
        double m_Eta;

        std::vector<double> m_Psi;       // Shape function values
        std::vector<double> m_PsiDKsi;   // Derivatives per ksi
        std::vector<double> m_PsiDEta;   // Derivatives per eta
    };

    //////////////////////////////////////////////////////////////////////////////////////////////
    ////  ILocalShapeFunctions1DLine
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Interface for local 1D line
    class ILocalShapeFunctions1DLine
    {
    public:
        virtual ~ILocalShapeFunctions1DLine() = default;

        explicit ILocalShapeFunctions1DLine(const LocalPoint1D & t_Point);

        virtual double Psi(size_t Index) const final;

    protected:
        double m_Ksi;
        std::vector<double> m_Psi;
    };

}   // namespace HygroThermFEM
