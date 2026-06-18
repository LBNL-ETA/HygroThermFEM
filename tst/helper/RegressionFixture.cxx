#include <utility>

#include "RegressionFixture.hxx"

namespace lbnl::errorest::test
{
    std::vector<Maker> & registry()
    {
        static std::vector<Maker> reg;
        return reg;
    }

    Registrar::Registrar(Maker maker)
    {
        registry().push_back(std::move(maker));
    }
}
