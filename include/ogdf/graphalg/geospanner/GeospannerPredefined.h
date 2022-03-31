/** \file
 * \brief Contains valid typedefs for geospanners
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

#include <ogdf/graphalg/geospanner/Geospanner.h>
#include <ogdf/graphalg/geospanner/SpannerDeltaGreedy.h>
#include <ogdf/graphalg/geospanner/SpannerPathGreedy.h>
#include <ogdf/graphalg/geospanner/SpannerYaoGraph.h>
#include <ogdf/graphalg/geospanner/geoutils/EuclidianAngleFinder.h>
#include <ogdf/graphalg/geospanner/geoutils/EuclidianMetric.h>
#include <ogdf/graphalg/geospanner/geoutils/HaversineMetric.h>
#include <ogdf/graphalg/geospanner/geoutils/SLCMetric.h>
#include <ogdf/graphalg/geospanner/geoutils/SphereAngleFinder.h>

namespace ogdf {

typedef SpannerDeltaGreedy<EuclidianMetric, EuclidianAngleFinder> SpannerDeltaGreedyEuclidian;
typedef SpannerPathGreedy<EuclidianMetric> SpannerPathGreedyEuclidian;
typedef SpannerYaoGraph<EuclidianMetric, EuclidianAngleFinder> SpannerYaoGraphEuclidian;

typedef SpannerDeltaGreedy<SLCMetric, SphereAngleFinder> SpannerDeltaGreedySphere;
typedef SpannerPathGreedy<SLCMetric> SpannerPathGreedySphere;
typedef SpannerYaoGraph<SLCMetric, SphereAngleFinder> SpannerYaoGraphSphere;

}  // namespace ogdf