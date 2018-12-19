#pragma once

#include <vector>
#include <memory>

#include "LocalShapeFunctions.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ////   QuadrilateralLinearLocal2D
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Represents quadrilateral linear element in local coordinate system (-1, 1)
    //!
    //! Note that class is done as singleton for simple reason that quadrilateral linear element
    //! in local coordinate system is always the same with no matter what coordinates of global
    //! element are.
    class QuadrilateralLinearLocal2D
    {
    public:
        static QuadrilateralLinearLocal2D & Instance();

        //! Local shape function
        double Psi(size_t IntegrationPointIndex,   //!< Index of integration point
                   size_t Index                    //!< Index of shape function
        );

        //! Differential of shape function over ksi (equivalent of x coordinate in local coordinate
        //! system)
        double PsiDKsi(size_t IntegrationPointIndex,   //!< Index of integration point
                       size_t Index                    //!< Index of shape function derivative
        );

        //! Differential of shape function over eta (equivalent of y coordinate in local coordinate
        //! system)
        double PsiDEta(size_t IntegrationPointIndex,   //!< Index of integration point
                       size_t Index                    //!< Index of shape function derivative
        );

        //! Vector all all shape functions in given integration point
		const std::vector< double > & Psi(
			size_t IntegrationPointIndex //!< Integration point index
		) const;

        //! Vector of all shape function derivatives over ksi (equivalent to x) in given integration
        //! point
		const std::vector< double > & PsiDKsi(
			size_t IntegrationPointIndex //!< Integration point index
		) const;

        //! Vector of all shape function derivatives over eta (equivalent to y) in given integration
        //! point
		const std::vector< double > & PsiDEta(
			size_t IntegrationPointIndex //!< Integration point index
		) const;

    private:
        QuadrilateralLinearLocal2D();

        //! Shape functions for every integration point (at local coordinate system)
        std::vector<std::unique_ptr<ILocalShapeFunctions2DQuadrilateral>> m_Ksi;

        //! \brief Class defines shape functions and its derivatives in local coordinate system
        //! \brief at given point.
        class QuadLinearLocalShapeFunctions2D : public ILocalShapeFunctions2DQuadrilateral
        {
        public:
            explicit QuadLinearLocalShapeFunctions2D(const LocalPoint2D & t_Point);
        };
    };

}   // namespace MoisThermFEM