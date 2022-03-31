/** \file
 * \brief Angle finder for points on a unit sphere for geospanners.
 * 
 * Projects points into azimuthal equidistant projection for angle calculation
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
 * @brief Objects to find deterministic angles from point to other points on spheres
 * Nodes must be associated with a GraphAttributes object with nodeGraphics enabled
 * 
 */
class SphereAngleFinder
{
public:
    SphereAngleFinder(const GraphAttributes &GA);

    /**
     * @brief Calculates an angle from a reference angle to the respective other point
     * 
     * First calculates the coordinates of each point in the azimuthal equidistant projection of the other point,
     * then uses the euclidian angle to the y-axis
     * 
     * Source on methodology used: https://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
     * 
     * Warning: latitudes must be between -p/2 and pi/2,longitudes between -pi and pi
     * 
     * @param v first point, given as ogdf::node
     * @param w second point, given as ogdf::node
     * @param distance 
     * @param eps 
     * @return the first entry is the (projected) angle of w from the perspective of v, the second 
     * the (projected) angle of v from the perspective of w. Both are between 0 and 2*Pi from (0,1) clockwise
     */
    std::pair<double, double> operator()(node v, node w, double distance);

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
    EpsilonTest m_eps;

    //stores sin and cos of latitude (eg. y of the respective node)
    NodeArray<std::pair<double, double>> preprocessed;
};

}  // namespace ogdf