#pragma once

#include <vector>
#include <ida/ida.h>                /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sundials/sundials_types.h>    /* definition of type realtype          */
#include <sunlinsol/sunlinsol_spgmr.h>  /* access to spgmr SUNLinearSolver      */
#include <sunmatrix/sunmatrix_band.h>   /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>   /* access to band SUNLinearSolver       */
#include "SolutionVariables.hxx"

namespace HygroThermFEM
{
    struct SingleDomain;
}

namespace Sundials
{
    // Put vector for now, but this is not going to work in long term because we need to switch
    // between different types of domains.
    std::vector<HygroThermFEM::SingleTimestepSolution>
      transient(HygroThermFEM::SingleDomain & domain,
                const std::vector<double> & previousTimestepValues,
                double t_DTime,
                size_t nTimesteps);
}   // namespace Sundials