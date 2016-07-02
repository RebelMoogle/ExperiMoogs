#pragma once

#include <Eigen\Core>
#include <boost/math/constants/constants.hpp>
// DoubleGyre

namespace FlowData {
    const static double pi = boost::math::constants::pi<double>(); // alternative, use M_PI. but since we're using boost anyway, this is more portable.

                                                                   // input: input position, output: output position, t: time
    static Eigen::Vector3d DoubleGyre3D(const Eigen::Vector3d& input, Eigen::Vector3d& output, const double t)
    {
        double x = input.x();
        double y = input.y();
        double z = input.z();

        output(0) = -(1.0 / 10.0)*std::sin((1.0 / 4.0)*pi*x*(std::sin((1.0 / 5.0)*pi*t)*x + 4.0 - 2.0 * std::sin((1.0 / 5.0)*pi*t)))*std::cos(pi*y)*pi;

        output(1) = (1.0 / 20.0)*pi*std::cos((1.0 / 4.0)*pi*x*(std::sin((1.0 / 5.0)*pi*t)*x + 4.0 - 2.0 * std::sin((1.0 / 5.0)*pi*t)))*(std::sin((1.0 / 5.0)*pi*t)*x + 2.0 - std::sin((1.0 / 5.0)*pi*t))*std::sin(pi*y);

        output(2) = (1.0 / 5.0)*z*(1.0 - z)*(z - (1.0 / 4.0)*std::sin((2.0 / 5.0)*pi*t) - 1.0 / 2.0);

        return output;
    }
}



