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

        //! \brief Creates and adds element into domain.
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        virtual void createElement(size_t index1,
                                   size_t index2,
                                   size_t index3,
                                   size_t index4,
                                   const std::string & materialName) override;

    protected:
        void postProcess(std::vector<double> & solution) override;
    };
}   // namespace HygroThermFEM
