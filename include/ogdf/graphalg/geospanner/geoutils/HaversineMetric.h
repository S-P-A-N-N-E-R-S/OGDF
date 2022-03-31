/** \file
 * \brief Metric for points on a unit sphere for geospanners.
 * 
 * Uses the haversine formular.
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

#pragma once

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/geometry.h>

namespace ogdf {

/**
 * @brief Calculates distance between two points by using the haversine formular
 * Nodes must be associated with a GraphAttributes object with nodeGraphics enabled
 */
class HaversineMetric
{
public:
    /**
     * @brief Construct a new Haversine Metric object
     * 
     * @param points the node array containing the coordinates where x is longitude and y latitude in radian
     */
    HaversineMetric(const GraphAttributes &graph);

    /**
     * @brief Calculates the distance between v and w on a unit sphere.
     * 
     * @param v 
     * @param w 
     * @return distance
     */
    inline double operator()(node v, node w) const
    {
        DPoint v_point{m_GA.x(v), m_GA.y(v)};
        DPoint w_point{m_GA.x(w), m_GA.y(w)};
        auto dif_point = w_point - v_point;
        auto lat_sin_dif = pow(sin(dif_point.m_y / 2), 2);
        double sum =
            lat_sin_dif + (1 - lat_sin_dif - pow(sin((v_point.m_y + w_point.m_y) / 2), 2)) *
                              pow(sin(dif_point.m_x / 2), 2);
        return 2 * asin(sqrt(sum));
    }

    /**
     * @brief Checks, if the GraphAttributes object can be used for this distance finder
     * 
     * @param GA The GraphAttributes Object to be checked
     * @param error If a precondition is hurt, this will be written onto error
     * @return 
     */
    static bool preconditionsOk(const GraphAttributes &GA, std::string &error);

private:
    //const NodeArray<DPoint> &points;
    const GraphAttributes &m_GA;
};

}  // namespace ogdf
