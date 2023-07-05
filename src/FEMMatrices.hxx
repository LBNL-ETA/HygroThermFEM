#pragma once

#include "SquareMatrix.hxx"
#include "SingleDomain.hxx"

namespace HygroThermFEM
{
    //! Forms left hand side matrix in steady state solution.
    //! @param domain Domain for which matrix is being formed
    SquareMatrix steadyStateLeftHandSide(SingleDomain & domain);

    //! Form right hand side vector in stead state solution.
    //! @param domain Domain for which vector is being formed
    [[nodiscard]] std::vector<double> steadyStateRightHandSide(SingleDomain & domain);

    //! \brief Forms mass, conductance and H (from boundary condition) matrices.
    //! @param domain Domain for which matrices are being formed
    //! @param t_DTime Time between two timestep for which calculates are being performed
    //! @param timestepIndex Timestep index used in case of variable timestep input boundary
    //! conditions.
    SquareMatrix transientM_K_H_Matrix(SingleDomain & domain, double t_DTime, size_t timestepIndex);

    //! \brief This function retrieves M*U+R vector (where U is state variable)
    //! @param domain Domain for which vector is being formed
    //! @param t_PreviousTimestepValues Solution from previous timestep
    //! @param t_DTime Time between two timestep for which calculates are being performed
    //! @param timestepIndex Timestep index used in case of variable timestep input boundary
    //! conditions.
    std::vector<double> transientMT_R_Vector(SingleDomain & domain,
                                             const std::vector<double> & t_PreviousTimestepValues,
                                             double t_DTime,
                                             size_t timestepIndex);

}   // namespace HygroThermFEM