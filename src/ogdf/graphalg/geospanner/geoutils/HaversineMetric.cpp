/** \file
 * \brief Implementation of haversine formular
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

#include <ogdf/graphalg/geospanner/geoutils/HaversineMetric.h>

namespace ogdf {

HaversineMetric::HaversineMetric(const GraphAttributes &graph)
    : m_GA(graph)
{
}

bool HaversineMetric::preconditionsOk(const GraphAttributes &GA, std::string &error)
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