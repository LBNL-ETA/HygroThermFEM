#pragma once

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint1D
    ////////////////////////////////////////////////////////////////////////////
    struct LocalPoint1D
    {
        explicit LocalPoint1D(double t_ksi);
        LocalPoint1D(const LocalPoint1D & t_LocalPoint);
        double ksi{0};
    };

    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint2D
    ////////////////////////////////////////////////////////////////////////////
    struct LocalPoint2D
    {
        LocalPoint2D(double t_ksi, double t_eta);
        LocalPoint2D(const LocalPoint2D & t_LocalPoint);
        double ksi{0};
        double eta{0};
    };
}   // namespace HygroThermFEM
