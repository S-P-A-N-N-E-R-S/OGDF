
/** \file
 * \brief Metric for points in 2D euclidian spaces for geospanners.
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
 * @brief Distance Finder of euclidian distance between two points.
 * Nodes must be associated with a GraphAttributes object with nodeGraphics enabled
 */
class EuclidianMetric
{
public:
    /**
     * @brief Construct a new Euclidian_Metric object
     * 
     * @param points the node array containing the coordinates
     */
    EuclidianMetric(const GraphAttributes &graph);

    /**
     * @brief Calculates euclidian distance between two nodes
     * 
     * @param v 
     * @param w 
     * @return the distance
     */
    inline double operator()(node v, node w) const
    {
        return sqrt(pow(m_GA.x(v) - m_GA.x(w), 2.0) + pow(m_GA.y(v) - m_GA.y(w), 2.0));
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
