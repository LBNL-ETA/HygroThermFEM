#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>
#include <sstream>

namespace MoisThermFEM
{

class ConvergenceException : public std::runtime_error
    {
    public:
        ConvergenceException();
        virtual const char* what() const throw();

    private:
    static std::ostringstream cnvt;
    };

}
