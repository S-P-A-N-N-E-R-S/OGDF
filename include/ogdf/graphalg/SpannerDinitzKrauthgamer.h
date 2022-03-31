/** \file
 * \brief Implementation of an r-vertex-fault-tolerant k-spanner algorithm from Dinitz
 * and Krauthgamer
 *
 * \author Julian Pascal Wittker
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

#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Algorithm for calculating fault-tolerant-spanners.
 *
 * Dinitz, M., & Krauthgamer, R. (2011, June).
 * Fault-tolerant spanners: better and simpler.
 * In Proceedings of the 30th annual ACM SIGACT-SIGOPS symposium
 * on Principles of distributed computing (pp. 169-178).
 *
 * Conditions for the graph:
 * - simple
 * - undirected
 *
 * The stretch mus satisfy \f$k\geq3\f$.
 *
 * The preconditions can be checked with SpannerDinitzKrauthgamer::preconditionsOk.
 *
 * Calculates a r-vertex-fault-tolerant k-spanner with \f$r\geq1\f$
 *
 * @ingroup ga-spanner <- TODO: ?
 */
template <typename TWeight, typename TSpanner>
class SpannerDinitzKrauthgamer : public SpannerModule<TWeight>
{
    static_assert(std::is_base_of<SpannerModule<TWeight>, TSpanner>::value,
                  "TSpanner must be a derived class of SpannerModule<TWeight> in "
                  "SpannerDinitzKrauthgamer<TWeight, TSpanner>.");

public:
    SpannerDinitzKrauthgamer()
        : m_r(0)
        , m_seed(-1)
    {
        setFaultSize(m_r);
        setIterScale(1);
    }

    void setFaultSize(int r)
    {
        m_r = r;
        m_r_Cubed = pow(m_r, 3);
        if (r == 1)
        {
            m_prob = 1.0 - 1.0 / 2.0;
        }
        else
        {
            m_prob = 1.0 - 1.0 / m_r;
        }
    }

    void setIterScale(int scale)
    {
        m_iter_scale = scale;
    }

    void setSeed(int seed)
    {
        m_seed = seed;
    }

    //! @copydoc ogdf::SpannerModule::preconditionsOk
    virtual bool preconditionsOk(const GraphAttributes &GA, double stretch,
                                 std::string &error) override
    {
        if (!isSimple(GA.constGraph()))
        {
            error = "The graph is not simple";
            return false;
        }
        if (stretch < 3.0)
        {
            error = "The stretch must be >= 3.0";
            return false;
        }
        if (GA.directed())
        {
            error = "The graph must be undirected";
            return false;
        }
        return true;
    }

private:
    const Graph *m_G;  //!< const reference to the original graph

    // the fault-tolerance parameter
    int m_r;
    // seed for random operations
    int m_seed;
    // the scale of the amount of iterations
    int m_iter_scale;

    // precalculated values
    int m_r_Cubed;   //!< r^3
    double m_n_Log;  //!< log(n)
    double m_prob;   //!< 1 - 1/r
    double m_iter;   //!< r^3 * log(n)

    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G = &GA.constGraph();
        m_n_Log = log(m_G->numberOfNodes());
        m_iter = m_r_Cubed * m_n_Log;

        if (m_seed == -1)
        {
            setSeed(time(NULL));
        }
        else
        {
            setSeed(m_seed);
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        // for o(r^3log(n)) iterations
        for (int i = 0; i < m_iter * m_iter_scale; i++)
        {
            assertTimeLeft();

            GraphCopySimple induced_graph;
            induced_graph.createEmpty(*m_G);
            List<node> selected_nodes;

            // remove vertices with probability 1-1/r from graph copy
            for (node n : m_G->nodes)
            {
                if (randomDouble(0, 1) > m_prob)
                {
                    selected_nodes.pushBack(n);
                }
            }
            inducedSubGraph(*m_G, selected_nodes.begin(), induced_graph);

            // H = H union k-Spanner(G\J)
            GraphCopySimple spanner_induced_graph(induced_graph);
            EdgeArray<bool> e_Arr(induced_graph);

            TSpanner induced_spanner;
            if (this->isTimelimitEnabled())
            {
                induced_spanner.setTimelimit(this->getTimeLeft());
            }
            GraphAttributes induced_ga(induced_graph, m_GA->attributes());
            induced_ga.directed() = false;
            for (edge e : induced_graph.edges)
            {
                induced_ga.doubleWeight(e) = m_GA->doubleWeight(induced_graph.original(e));
            }
            auto res = induced_spanner.call(induced_ga, m_stretch, spanner_induced_graph, e_Arr);

            // adapt spanner's edge array
            for (edge e : induced_graph.edges)
            {
                if (e_Arr[e])
                {
                    node a = e->nodes()[0];
                    node b = e->nodes()[1];
                    if (m_spanner->searchEdge(m_spanner->copy(induced_graph.original(a)),
                                              m_spanner->copy(induced_graph.original(b))) ==
                        nullptr)
                    {
                        m_spanner->newEdge(induced_graph.original(e));
                        (*m_inSpanner)[induced_graph.original(e)] = true;
                    }
                }
            }
        }

        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    using SpannerModule<TWeight>::getWeight;
    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_stretch;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};

}  // namespace ogdf
