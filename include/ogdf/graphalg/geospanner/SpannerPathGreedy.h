/** \file
 * \brief Implementation of basic path greedy geospanner
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
#include <ogdf/basic/NodeSet.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * @brief Path greedy algorithm for geospanner. A description can be found
 * in "Geometric Spanner Networks" by Narasimhan and Smid.
 * 
 * @tparam Metric 
 */
template <class Metric>
class SpannerPathGreedy : public Geospanner<Metric>
{
    using TDistance = TMetricDistance<Metric>;

    EdgeArray<TDistance> m_spannerWeights;
    EpsilonTest m_eps;

public:
    //! @copydoc ogdf::SpannerModule::preconditionsOk
    virtual bool preconditionsOk(const GraphAttributes &GA, double stretch,
                                 std::string &error) override
    {
        if (!GA.has(GraphAttributes::nodeGraphics))
        {
            error = "Node coordinates must be provided via graph attributes";
            return false;
        }
        if (stretch < 1.0)
        {
            error = "The stretch must be >= 1.0";
            return false;
        }
        if (!Metric::preconditionsOk(GA, error))
        {
            return false;
        }
        return true;
    }

    //! @copydoc ogdf::Geospanner::weights
    const EdgeArray<TDistance> &weights()
    {
        return m_spannerWeights;
    }

private:
    using SpannerModule<TDistance>::assertTimeLeft;
    using SpannerModule<TDistance>::m_GA;
    using SpannerModule<TDistance>::m_stretch;
    using SpannerModule<TDistance>::m_spanner;
    using typename Geospanner<Metric>::edgeEntry;
    using typename Geospanner<Metric>::edgeEntryComparator;

    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TDistance>::init(GA, stretch, spanner, inSpanner);
        m_spannerWeights.init(spanner);
    }

    virtual typename SpannerModule<TDistance>::ReturnType execute() override
    {
        Metric distance(*m_GA);

        const Graph &og_graph = m_GA->constGraph();

        // Create an array for all potential edges and store them in it
        Array<edgeEntry> edges(og_graph.numberOfNodes() * (og_graph.numberOfNodes() - 1) / 2);

        int i = 0;
        for (auto u = og_graph.nodes.begin(); u != og_graph.nodes.end(); u++)
        {
            auto v = u;
            v++;
            for (; v != og_graph.nodes.end(); v++)
            {
                edges[i++] = edgeEntry{*u, *v, distance(*u, *v)};
            }
            assertTimeLeft();
        }

        // Sort potential edges
        edges.quicksort(edgeEntryComparator());

        assertTimeLeft();

        // Use greedy approach to decide which edge is inserted
        for (const edgeEntry &e : edges)
        {
            node u = m_spanner->copy(e.u);
            node v = m_spanner->copy(e.v);
            double maxDistance = m_stretch * e.distance;
            double currentSpannerDistance = distanceInSpanner(u, v, maxDistance);
            if (m_eps.greater(currentSpannerDistance, maxDistance))
            {
                edge newEdge = m_spanner->newEdge(u, v);
                m_spannerWeights[newEdge] = e.distance;
            }
            assertTimeLeft();
        }

        return SpannerModule<TDistance>::ReturnType::Feasible;
    }

    double distanceInSpanner(node s, node t, double maxLookupDist)
    {
        NodeArray<TDistance> distances;
        NodeArray<edge> predecessor;
        Dijkstra<TDistance> dijkstra;
        dijkstra.callBound(*m_spanner, m_spannerWeights, s, predecessor, distances,
                           false,  // directed
                           false,  // arcs reversed
                           t, maxLookupDist);
        return distances[t];
    }
};

}  // namespace ogdf
