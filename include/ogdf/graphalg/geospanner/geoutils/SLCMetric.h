/** \file
 * \brief Metric for points on a unit sphere for geospanners.
 * 
 * Uses the spherical law of cosines. Generally faster than haversine formular, but more prone
 * to numerical errors if distances are very low.
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
 * @brief Calculates distance between two points by using the spherical law of cosines
 * Nodes must be associated with a GraphAttributes object with nodeGraphics enabled
 */
class SLCMetric
{
public:
    /**
     * @brief Construct a new Haversine Metric object
     * 
     * @param points the node array containing the coordinates where x is longitude and y latitude in radian
     */
    SLCMetric(const GraphAttributes &graph);

    /**
     * @brief Calculates the distance between v and w on a unit sphere.
     * 
     * @param v 
     * @param w 
     * @return distance
     */
    inline double operator()(node v, node w) const
    {
        double v_lon = m_GA.x(v);
        double w_lon = m_GA.x(w);
        auto v_pre = preprocessed[v];
        auto w_pre = preprocessed[w];
        return acos(v_pre.first * w_pre.first + v_pre.second * w_pre.second * cos(w_lon - v_lon));
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
    const GraphAttributes &m_GA;

    //stores sin and cos of latitude (eg. y of the respective node)
    NodeArray<std::pair<double, double>> preprocessed;
};

}  // namespace ogdf