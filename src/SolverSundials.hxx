#pragma once

#include <vector>
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
                size_t nTimesteps = 0u);
}   // namespace Sundials