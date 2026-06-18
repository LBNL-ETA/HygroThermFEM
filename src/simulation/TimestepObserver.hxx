#pragma once

#include "Common.hxx"

namespace Timesteps {

    class TimestepObserver {
    public:
        virtual ~TimestepObserver() = default;

        virtual void levelChanged(unsigned divisionLevel, unsigned timestepNumber) = 0;
    };

}
