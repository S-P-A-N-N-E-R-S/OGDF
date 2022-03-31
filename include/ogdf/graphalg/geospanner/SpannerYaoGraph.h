/** \file
 * \brief Implementation of basic yao graph algorithm by Yao (1982)
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
#include <ogdf/graphalg/geospanner/Geospanner.h>
#include <cmath>
#include <limits>
#include <vector>

namespace ogdf {

/**
 * Building a geometric spanner via choosing nearest neighbour in cones
 * 
 * Described by Yao in "On constructing minimum spanning trees in k-dimensional space
 * and related problems" (1982)
 * 
 * Conditions: Coordinates for vertices.
 * 
 */
template <class Metric, class AngleFinder>
class SpannerYaoGraph : public Geospanner<Metric>
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
        return true;
    }

    //! @copydoc ogdf::Geospanner::weights
    virtual const EdgeArray<TDistance> &weights() override
    {
        return m_weights;
    }

    /**
     * @brief Stores the information to which node in a cone an edge exists
     * 
     */
    struct ConeEntry {
        /**
         * @brief Distance to the reference node in this cone
         * 
         */
        TDistance distance;

        /**
         * @brief Reference node in this cone
         * 
         */
        node target;

        /**
         * @brief Index of the cone of the reference node this node would be in
         * 
         */
        int target_cone_index;
    };

    const NodeArray<std::vector<ConeEntry>> &node_cones()
    {
        return m_node_cones;
    }

private:
    using SpannerModule<TDistance>::assertTimeLeft;
    using SpannerModule<TDistance>::m_GA;
    using SpannerModule<TDistance>::m_stretch;
    using SpannerModule<TDistance>::m_spanner;

    NodeArray<std::vector<ConeEntry>> m_node_cones;
    EdgeArray<TDistance> m_weights;
    EpsilonTest m_eps;

    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TDistance>::init(GA, stretch, spanner, inSpanner);
        m_node_cones.init(spanner);
        m_weights.init(spanner);
    }

    virtual typename SpannerModule<TDistance>::ReturnType execute() override
    {
        // Reserve memory for cones and initialize cones as empty

        int k = ceil(2 * M_PI / (1 - 1 / m_stretch));
        const Graph &og_graph = m_GA->constGraph();
        for (auto &node_cone : m_node_cones)
        {
            node_cone = std::vector<ConeEntry>(
                k, ConeEntry{std::numeric_limits<TDistance>::max(), nullptr, -1});
        }

        assertTimeLeft();

        Metric distance(*m_GA);
        AngleFinder angle(*m_GA);

        assertTimeLeft();

        for (auto v = og_graph.nodes.begin(); v != og_graph.nodes.end(); v++)
        {
            auto v_spanner = m_spanner->copy(*v);
            auto &cones_v = m_node_cones[v_spanner];

            // we iterate over all node pairs once
            auto w = v;
            w++;
            for (; w != og_graph.nodes.end(); w++)
            {
                TDistance distance_vw = distance(*v, *w);

                auto angles = angle(*v, *w, distance_vw);

                double angle_v = angles.first;
                int cone_index_v = angle_v * k / (2 * M_PI);
                if (cone_index_v == k)  // can happen due to numerical errors
                {
                    cone_index_v = 0;
                }
                double angle_w = angles.second;
                int cone_index_w = angle_w * k / (2 * M_PI);
                if (cone_index_w == k)  // can happen due to numerical errors
                {
                    cone_index_w = 0;
                }

                auto w_spanner = m_spanner->copy(*w);

                if (m_eps.less(distance_vw, cones_v[cone_index_v].distance))
                {
                    cones_v[cone_index_v] = ConeEntry{distance_vw, w_spanner, cone_index_w};
                }

                auto &cones_w = m_node_cones[w_spanner];
                if (m_eps.less(distance_vw, cones_w[cone_index_w].distance))
                {
                    cones_w[cone_index_w] = ConeEntry{distance_vw, v_spanner, cone_index_v};
                }
            }
            assertTimeLeft();
        }

        for (auto v : m_spanner->nodes)
        {
            auto node_cone = m_node_cones[v];
            for (const auto &t : node_cone)
            {
                if (t.target != nullptr)
                {
                    if (t.target->index() < v->index() &&
                        m_node_cones[t.target][t.target_cone_index].target == v)
                    {
                        // If target has a smaller index and v is also its reference node in a cone, then the edge is already build by target
                        continue;
                    }
                    auto e = m_spanner->newEdge(v, t.target);
                    m_weights[e] = t.distance;
                }
            }
        }

        return SpannerModule<TDistance>::ReturnType::Feasible;
    }
};

}  // namespace ogdf