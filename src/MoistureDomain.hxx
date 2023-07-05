#pragma once

#include "SingleDomain.hxx"

namespace HygroThermFEM
{
    //! \brief Domain class for solving humidity distribution.
    struct MoistureDomain : public SingleDomain
    {
    public:
        //! Simple constructor
        MoistureDomain(bool automaticUpdatePreviousTimestep = true);

    protected:
        void postProcess(std::vector<double> & solution) override;
    };
}   // namespace HygroThermFEM
