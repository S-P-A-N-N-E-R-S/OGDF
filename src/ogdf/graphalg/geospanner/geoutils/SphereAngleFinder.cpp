/** \file
 * \brief Implementation of spherical angle finder
 *
 * \author Levin Nemesch
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */


#include <ogdf/graphalg/geospanner/geoutils/SphereAngleFinder.h>

namespace ogdf {

SphereAngleFinder::SphereAngleFinder(const GraphAttributes &GA)
    : m_GA(GA)
    , preprocessed(m_GA.constGraph())
{
    // Preprocesses trigonometric functions for latitudes
    for (auto v : m_GA.constGraph().nodes)
    {
        double lat = m_GA.y(v);
        preprocessed[v] = {sin(lat), cos(lat)};
    }
}

std::pair<double, double> SphereAngleFinder::operator()(node v, node w, double distance)
{
    auto lon_v = m_GA.x(v);
    //auto v_pre = preprocessed[v];
    double sin_lat_v, cos_lat_v;
    std::tie(sin_lat_v, cos_lat_v) = preprocessed[v];

    auto lon_w = m_GA.x(w);
    //auto w_pre = preprocessed[w];
    double sin_lat_w, cos_lat_w;
    std::tie(sin_lat_w, cos_lat_w) = preprocessed[w];

    if (m_eps.equal(distance, M_PI))
    {
        return {0.0, 0.0};
    }

    double k = distance / sin(distance);
    double sin_w_minus_v = sin(lon_w - lon_v);
    double cos_w_minus_v = cos(lon_w - lon_v);

    // angle of w from perspective of v
    double x_v = k * cos_lat_w * sin_w_minus_v;
    double y_v = k * (cos_lat_v * sin_lat_w - sin_lat_v * cos_lat_w * cos_w_minus_v);

    // switch x and y in atan2 so we measure angle to (0,1)
    double angle_v = m_eps.equal(x_v, 0.0) && m_eps.equal(y_v, 0.0) ? 0.0 : atan2(x_v, y_v);

    // angle of v from perspective of w
    double x_w = k * cos_lat_v * (-sin_w_minus_v);
    double y_w = k * (cos_lat_w * sin_lat_v - sin_lat_w * cos_lat_v * cos_w_minus_v);

    // switch x and y in atan2 so we measure angle to (0,1)
    double angle_w = m_eps.equal(x_w, 0.0) && m_eps.equal(y_w, 0.0) ? 0.0 : atan2(x_w, y_w);

    return {angle_v < 0 ? M_PI * 2 + angle_v : angle_v, angle_w < 0 ? M_PI * 2 + angle_w : angle_w};
}

bool SphereAngleFinder::preconditionsOk(const GraphAttributes &GA, std::string &error)
{
    if (!GA.has(GraphAttributes::nodeGraphics))
    {
        error = "Node Coordinates must be active in GraphAttributes";
        return false;
    }
    const EpsilonTest epsilon;
    for (const auto n : GA.constGraph().nodes)
    {
        if (epsilon.less(M_PI, fabs(GA.x(n))))
        {
            error = "All longitudes must be between -pi and pi";
            return false;
        }
        if (epsilon.less(M_PI_2, fabs(GA.y(n))))
        {
            error = "All latitudes must be between -pi/2 and pi/2";
            return false;
        }
    }
    return true;
}

}  // namespace ogdf