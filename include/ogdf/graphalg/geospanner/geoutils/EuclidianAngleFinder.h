/** \file
 * \brief Angle finder for points in 2D euclidian spaces for geospanners.
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
#include <cmath>
#include <utility>

namespace ogdf {

/**
 * @brief Angle Finder for two euclidian points.
 * Nodes must be associated with a GraphAttributes object with nodeGraphics enabled
 */
class EuclidianAngleFinder
{
public:
    /**
     * @brief Construct a new Euclidian_AngleFinder object
     * 
     * @param points the node array containing the coordinates
     */
    EuclidianAngleFinder(const GraphAttributes &GA);

    /**
     * @brief Calculates the respective angle of (0,1) to w-v (for angle of w to v) and v-w
     * 
     * @param v 
     * @param w 
     * @param distance 
     * @return the angles as radians between 0 and 2*Pi
     */
    inline std::pair<double, double> operator()(node v, node w, double distance)
    {
        auto x_v = m_GA.x(v);
        auto y_v = m_GA.y(v);

        auto x_w = m_GA.x(w);
        auto y_w = m_GA.y(w);

        // switch x and y in atan2 so we measure angle to (0,1)
        double angle_v = atan2(x_w - x_v, y_w - y_v);
        double angle_w = atan2(x_v - x_w, y_v - y_w);
        return {angle_v < 0 ? M_PI * 2 + angle_v : angle_v,
                angle_w < 0 ? M_PI * 2 + angle_w : angle_w};
    }

    /**
     * @brief Checks, if the GraphAttributes object can be used for this angle finder
     * 
     * @param GA The GraphAttributes Object to be checked
     * @param error If a precondition is hurt, this will be written onto error
     * @return 
     */
    static bool preconditionsOk(const GraphAttributes &GA, std::string &error);

private:
    const GraphAttributes &m_GA;
};

}  // namespace ogdf