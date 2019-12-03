#pragma once

#include "Common.hxx"

template <typename T>
class TimestepObserver {
public:
    virtual ~TimestepObserver() = default;
    virtual void fieldChanged(const Timestep::Level& timestepLevel, unsigned timestepNumber) = 0;
};
