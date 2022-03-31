/** \file
 * \brief Base class for geospanner algorithms.
 * 
 * Every geospanner derives from SpannerModule and from Geospanner class.
 * 
 * The function weight provides access to the weights of the edges in the resulting spanner.
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>
#include <ogdf/graphalg/SpannerModule.h>
#include <utility>

namespace ogdf {

/**
 * @brief Specifies return type of the distance measurement of a Metric
 * 
 * @tparam Metric 
 */
template <class Metric>
using TMetricDistance =
    decltype(std::declval<Metric>().operator()(std::declval<node>(), std::declval<node>()));

/**
 * @brief Base class for geospanners. Also contains functions to check correctness and actual stretch factor
 * 
 * Different to other spanners, a geospanners only precondition are node coordinates and a stretch >= 1.
 * The inSpanner EdgeArray given to call is meaningless and completely ignored by geospanner algorithms
 * Furthermore, the option to access the edge weights of the spanner is provided by the Geospanner class.
 * SpannerModules isMultiplicativeSpanner wont work, isGeospanner provides this functionality.
 * Everything else behaves like defined in SpannerModule.h
 * 
 * @tparam Metric
 */
template <class Metric>
class Geospanner : public SpannerModule<TMetricDistance<Metric>>
{
public:
    using TDistance = TMetricDistance<Metric>;

    /**
     * @brief Returns the weights of the edges in the final spanner
     * 
     * @return const EdgeArray<TDistance>& 
     */
    virtual const EdgeArray<TDistance> &weights() = 0;

    /**
     * @brief Checks whether a graph is a geospanner to a given stretch
     * 
     * @tparam Metric 
     * @param GA Original GraphAttributes with the coordinates
     * @param spanner The GraphCopy with the spanner edges. Nodes need to correspond GA
     * @param weights Weights of the spanner edges
     * @param stretch stretch factor to test against
     * @return true if actual stretch factor is smaller or equal than stretch, false if not
     */
    static bool isGeospanner(const GraphAttributes &GA, const GraphCopySimple &spanner,
                             const EdgeArray<TDistance> &weights, double stretch)
    {
#ifdef OGDF_DEBUG
        std::string error;
        OGDF_ASSERT(Metric::preconditionsOk(GA, error));
        OGDF_ASSERT(&(spanner.original()) == &GA.constGraph());
#endif

        NodeArray<TDistance> spannerDistances(spanner);

        Metric distance(GA);
        EpsilonTest epsilon;

        for (auto v = spanner.nodes.begin(); v != spanner.nodes.end(); v++)
        {
            dijkstra_SPSS(*v, spanner, spannerDistances, weights);

            auto w = v;
            w++;

            for (; w != spanner.nodes.end(); w++)
            {
                if (epsilon.greater(spannerDistances[*w],
                                    stretch * distance(spanner.original(*v), spanner.original(*w))))
                {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * @brief Checks actual stretch factor of a spanner
     * 
     * @tparam Metric 
     * @param GA Original GraphAttributes with the coordinates
     * @param spanner The GraphCopy with the spanner edges. Nodes need to correspond GA
     * @param weights Weights of the spanner edges
     * @return actual stretch factor
     */
    static TDistance maxStretch(const GraphAttributes &GA, const GraphCopySimple &spanner,
                                const EdgeArray<TDistance> &weights)
    {
#ifdef OGDF_DEBUG
        std::string error;
        OGDF_ASSERT(Metric::preconditionsOk(GA, error));
        OGDF_ASSERT(&(spanner.original()) == &GA.constGraph());
#endif

        NodeArray<TDistance> spannerDistances(spanner);

        Metric distance(GA);
        EpsilonTest epsilon;
        TDistance max_stretch{1.0};

        for (auto v = spanner.nodes.begin(); v != spanner.nodes.end(); v++)
        {
            dijkstra_SPSS(*v, spanner, spannerDistances, weights);

            auto w = v;
            w++;

            for (; w != spanner.nodes.end(); w++)
            {
                TDistance test_val{spannerDistances[*w] /
                                   distance(spanner.original(*v), spanner.original(*w))};
                if (epsilon.greater(test_val, max_stretch))
                {
                    max_stretch = test_val;
                }
            }
        }

        return max_stretch;
    }

protected:
    /**
     * @brief Struct to store potential edges of the graph
     * 
     */
    struct edgeEntry {
        node u;
        node v;
        TDistance distance;
    };

    /**
     * @brief Compares two edgeEntry object by their respective distance
     * 
     */
    struct edgeEntryComparator {
        const EpsilonTest m_eps;
        bool less(edgeEntry a, edgeEntry b) const
        {
            return m_eps.less(a.distance, b.distance);
        }
    };
};

}  // namespace ogdf