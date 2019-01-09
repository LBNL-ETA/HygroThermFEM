#include "ConvergenceException.hxx"

namespace MoisThermFEM
{
    ConvergenceException::ConvergenceException() : std::runtime_error("Program failed to converge")
    {}

    const char * ConvergenceException::what() const throw()
    {
        cnvt.str( "" );

        cnvt << runtime_error::what();

        return cnvt.str().c_str();
    }
}   // namespace MoisThermFEM