/** \file
 * \brief Implementation of basic delta greedy algorithm by Abu-Affash, Bar-On and Carmi (2022)
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

/**
 * @file SpannerDeltaGreedy.h
 * @author Levin Nemesch (lnemesch@uos.de)
 * 
 * @date 2021-12-09
 */

#include <ogdf/basic/NodeSet.h>
#include <ogdf/graphalg/SpannerModule.h>
#include <cmath>
#include <limits>
#include <vector>

namespace ogdf {

/**
 * @brief Class for the delta greedy algorithm described by 
 * Abu-Affash, Bar-On and Carmi in "δ-Greedy t-spanner" (2022)
 * 
 * @tparam Metric 
 * @tparam AngleFinder 
 */
template <class Metric, class AngleFinder>
class SpannerDeltaGreedy : public Geospanner<Metric>
{
public:
    using TDistance = TMetricDistance<Metric>;

    //! @copydoc ogdf::SpannerModule::preconditionsOk
    virtual bool preconditionsOk(const GraphAttributes &GA, double stretch,
                                 std::string &error) override
    {
        if (!GA.has(GraphAttributes::nodeGraphics))
        {
            error = "Node coordinates must be provided via graph attributes";
            return false;
        }
        if (stretch <= 1.0)
        {
            error = "The stretch must be > 1.0";
            return false;
        }
        if (!Metric::preconditionsOk(GA, error))
        {
            return false;
        }
        if (!AngleFinder::preconditionsOk(GA, error))
        {
            return false;
        }
        if (m_eps.geq(m_delta, stretch))
        {
            error = "delta must be < t";
            return false;
        }
        return true;
    }

    //! @copydoc ogdf::Geospanner::weights
    virtual const EdgeArray<TDistance> &weights() override
    {
        return m_weights;
    }

    /**
     * @brief Sets delta. 
     * If delta is <1, default delta is chosen as =sqrt(stretch)
     * If delta is >stretch delta is chosen as =stretch
     * 
     * @param delta 
     */
    void setDelta(double delta)
    {
        m_delta = delta;
    }

private:
    using SpannerModule<TDistance>::assertTimeLeft;
    using SpannerModule<TDistance>::m_GA;
    using SpannerModule<TDistance>::m_stretch;
    using SpannerModule<TDistance>::m_spanner;
    using typename Geospanner<Metric>::edgeEntry;
    using typename Geospanner<Metric>::edgeEntryComparator;

    EdgeArray<TDistance> m_weights;
    double m_delta = 0;
    EpsilonTest m_eps;

    /**
     * @brief A cone is defined by its reference angle (on a line to another point) and its width
     * 
     */
    struct ConeEntryDelta {
        /**
         * @brief Reference angle in this cone
         * 
         */
        double reference_angle;

        /**
         * @brief half the width of the cone, eg the angle from the reference angle to an edge
         * 
         */
        double width_half;
    };

    NodeArray<std::vector<ConeEntryDelta>> m_node_cones;

    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TDistance>::init(GA, stretch, spanner, inSpanner);
        m_weights.init(spanner);
        m_node_cones.init(spanner);
        for (auto &node_cone : m_node_cones)
        {
            node_cone = std::vector<ConeEntryDelta>();
        }
    }

    virtual typename SpannerModule<TDistance>::ReturnType execute() override
    {
        // Delta always needs to be valid so we choose a valid delta if invalid
        double delta = m_eps.less(m_delta, 1.0) ? sqrt(m_stretch) : m_delta;

        Metric distance(*m_GA);
        AngleFinder angle(*m_GA);

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

        // sort edges
        edges.quicksort(edgeEntryComparator());

        assertTimeLeft();

        for (const edgeEntry &e : edges)
        {
            double angle_v, angle_w;
            std::tie(angle_v, angle_w) = angle(e.u, e.v, e.distance);
            node v = m_spanner->copy(e.u);
            node w = m_spanner->copy(e.v);
            if (inCone(v, angle_v) || inCone(w, angle_w))
            {
                continue;
            }

            TDistance maxDistance = m_stretch * e.distance;
            TDistance currentSpannerDistance = distanceInSpanner(v, w, maxDistance);

            TDistance d = currentSpannerDistance / e.distance;

            if (m_eps.greater(d, delta))
            {
                auto newEdge = m_spanner->newEdge(v, w);
                m_weights[newEdge] = e.distance;
                d = 1;
            }

            double width = M_PI_4 - asin(d / (sqrt(2) * m_stretch));

            m_node_cones[v].push_back(ConeEntryDelta{angle_v, width});
            m_node_cones[w].push_back(ConeEntryDelta{angle_w, width});

            assertTimeLeft();
        }

        return SpannerModule<TDistance>::ReturnType::Feasible;
    }

    /**
     * @brief Checks if an angle to a point is in one of its cones
     * 
     * @param v 
     * @param angle angle to check
     * @return true if angle is in an arbitrary cone, false otherwise
     */
    bool inCone(node v, double angle)
    {
        // Linear search faster for most realistically sized node cones
        for (const auto &nodeCone : m_node_cones[v])
        {
            double ref = nodeCone.reference_angle;
            double width = nodeCone.width_half;

            // Checks if w lays in the cone without taking overlap under 0 or above 2*Pi into account
            if (m_eps.geq(angle, ref - width) && m_eps.leq(angle, ref + width))
            {
                return true;
            }

            // If the cone overlaps into negative territory, checks if w is in an equivalent cone near 2*Pi
            if (m_eps.less(ref - width, 0.0) && m_eps.geq(angle, 2 * M_PI + ref - width))
            {
                return true;
            }

            // If the cone overlaps above 2*Pi, checks if w is in an equivalent cone near 0
            if (m_eps.greater(ref + width, 2 * M_PI) && m_eps.leq(angle, ref + width - 2 * M_PI))
            {
                return true;
            }
        }

        return false;
    }

    TDistance distanceInSpanner(node s, node t, TDistance maxLookupDist)
    {
        NodeArray<TDistance> distances;
        NodeArray<edge> predecessor;
        Dijkstra<TDistance> dijkstra;
        dijkstra.callBound(*m_spanner, m_weights, s, predecessor, distances,
                           false,  // directed
                           false,  // arcs reversed
                           t, maxLookupDist);
        return distances[t];
    }
};

}  // namespace ogdf