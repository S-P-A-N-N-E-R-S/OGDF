/** \file
 * \brief Implementation of a (2k-1)-spanner algorithm
 * from Thorup and Zwick
 *
 * \author Tim Hartmann, Leon Nienhüser
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

#include <ogdf/basic/Queue.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/SpannerModule.h>
#include <algorithm>

namespace ogdf {

/**
 * Algorithm for calculating (2k-1)-spanners.
 *
 * M. Thorup and U. Zwick.
 * Approximate distance oracles (2005).
 * J. ACM 52, 1 (January 2005), 1–24.
 * DOI:https://doi.org/10.1145/1044731.1044732
 *
 * Conditions for the graph:
 * - simple
 * - undirected
 * - weighted
 *
 * The stretch must satisfy \f$k\geq1\f$.
 *
 * The preconditions can be checked with SpannerThorupZwick::preconditionsOk.
 *
 * Calculates a (2k-1)-spanner with \f$k\geq1\f$
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class SpannerThorupZwick : public SpannerModule<TWeight>
{
public:
    SpannerThorupZwick()
        : m_seed(-1)
    {
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
        if (GA.directed())
        {
            error = "The graph must be undirected";
            return false;
        }
        int intStretch = static_cast<int>(stretch);
        if (intStretch < 1)
        {
            error = "The stretch must be >= 1.0";
            return false;
        }
        if (intStretch % 2 == 0)
        {
            error = "The stretch is required to be an integer, not " + to_string(stretch);
            return false;
        }
        double integralPart;
        if (std::modf((stretch + 1) / 2, &integralPart) != 0.0)
        {
            error = "The stretch is required to be an odd integer, not " + to_string(stretch);
            return false;
        }
        return true;
    }

private:
    // original graph
    GraphCopySimple m_G;
    // the parameter k derived from the stretch
    int m_k;
    // seed for random operations
    int m_seed;

    EdgeArray<TWeight> m_weights;

    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_k = (stretch + 1) / 2;
        m_G.init(GA.constGraph());

        m_weights.init(m_G);
        for (const edge e : m_G.edges)
        {
            m_weights[e] = SpannerModule<TWeight>::getWeight(GA, m_G.original(e));
        }

        if (m_seed == -1)
        {
            setSeed(time(NULL));
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        assertTimeLeft();

        std::vector<std::vector<node>> a(m_k + 1);
        a[0].resize(m_G.numberOfNodes());
        int count = 0;
        for (const node v : m_G.nodes)
        {
            a[0][count++] = v;
        }

        double exp = -1.0 / m_k;
        double probability = pow(m_G.numberOfNodes(), exp);

        for (int i = 1; i <= m_k - 1; i++)
        {
            for (const node n : a[i - 1])
            {
                double randomD = randomDouble(0.0, 1.0);
                if (randomD <= probability)
                {
                    a[i].push_back(n);
                }
            }
        }

        std::vector<NodeArray<TWeight>> delta;
        // init m_delta (distances of nodes to each of the witnesses)
        delta.resize(m_k + 1);
        delta[m_k].init(m_G, std::numeric_limits<TWeight>::max());

        NodeArray<std::vector<node>> clusters(m_G);

        // do k SSSP calls
        for (int i = m_k - 1; i >= 0; i--)
        {
            assertTimeLeft();

            node source = m_G.newNode();

            for (const node n : a[i])
            {
                edge e = m_G.newEdge(source, n);
                m_weights[e] = 0.0;
            }

            // compute SSSP for newly inserted source vertex
            NodeArray<edge> preds;
            NodeArray<TWeight> distResult;

            // call dijkstra to get shortest distances
            Dijkstra<TWeight>().call(m_G, m_weights, source, preds, distResult, false);

            delta[i] = distResult;

            m_G.delNode(source);

            for (const node w : a[i])
            {
                if (std::find(a[i + 1].begin(), a[i + 1].end(), w) == a[i + 1].end())
                {
                    modDijkstra(w, clusters[w], delta, i);
                }
            }
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
	 * Executes a modified dijkstra search to form clusters. The modification to the standard
	 * dijkstra algorithm is that we relax an edge (u,v) only if the new distance found for
	 * the vertex v is smaller then the \p delta value of v of the next \p iteration.
	 *
	 * @param source The start node of the search
	 * @param cluster Contains all the nodes visited by the search
	 * @param delta The precalculated distances to the samples
	 * @param iteration The current iteration
	 */
    void modDijkstra(const node source, std::vector<node> &cluster,
                     const std::vector<NodeArray<TWeight>> &delta, const int iteration)
    {
        assertTimeLeft();

        // run modified dijkstra
        PrioritizedMapQueue<node, TWeight, std::less<TWeight>> pq(m_G);
        NodeArray<TWeight> distanceReturn(m_G, std::numeric_limits<TWeight>::max());

        for (const node v : m_G.nodes)
        {
            pq.push(v, distanceReturn[v]);
        }

        distanceReturn[source] = 0.0;
        pq.decrease(source, distanceReturn[source]);
        while (!pq.empty())
        {
            node u = pq.topElement();
            pq.pop();
            if (distanceReturn[u] == std::numeric_limits<TWeight>::max())
            {
                break;
            }
            cluster.push_back(u);
            for (adjEntry adj : u->adjEntries)
            {
                edge e = adj->theEdge();
                node v = e->opposite(u);

                // check conditions for relaxation
                auto weight = SpannerModule<TWeight>::getWeight(*m_GA, m_G.original(e));

                if (distanceReturn[v] > distanceReturn[u] + weight &&
                    delta[iteration + 1][v] > distanceReturn[u] + weight)
                {
                    distanceReturn[v] = distanceReturn[u] + weight;
                    pq.decrease(v, distanceReturn[v]);
                    if (m_spanner->searchEdge(m_spanner->copy(m_G.original(u)),
                                              m_spanner->copy(m_G.original(v))) == nullptr)
                    {
                        m_spanner->newEdge(m_G.original(e));
                        (*m_inSpanner)[m_G.original(e)] = true;
                    }
                }
            }
        }
    }

    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};

}  // namespace ogdf
