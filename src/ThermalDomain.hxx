#pragma once

#include "SingleDomain.hxx"

namespace HygroThermFEM
{
    //! \brief Domain class for solving temperature solution.
    struct ThermalDomain : public SingleDomain
    {
        //! Simple constructor
        ThermalDomain(bool automaticUpdatePreviousTimestep = true);

    protected:
        void postProcess(std::vector<double> & solution) override;
    };
}   // namespace HygroThermFEM
