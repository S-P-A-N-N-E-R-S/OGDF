/** \file
 * \brief Implementation of an r-vertex-fault-tolerant (2k-1)-spanner algorithm
 * from Chechik and Langberg
 *
 * \author Tim Hartmann, Leon Nienhüser, Julian Pascal Wittker
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
 * Algorithm for calculating fault-tolerant-(2k-1)-spanners.
 *
 * Chechik, S., Langberg, M., Peleg, D., & Roditty, L. (2010).
 * Fault tolerant spanners for general graphs.
 * SIAM Journal on Computing, 39(7), 3403-3423.
 *
 * Conditions for the graph:
 * - simple
 * - undirected
 * - unweighted
 *
 * The stretch mus satisfy \f$k\geq1\f$.
 *
 * The preconditions can be checked with SpannerChechikLangbergVertexSimple::preconditionsOk.
 *
 * Calculates a r-edge-fault-tolerant (2k-1)-spanner with \f$r\geq1\f$
 * with inspiration from the Thorup Zwick Spanner
 *
 * @ingroup ga-spanner <- TODO: ?
 */
template <typename TWeight>
class SpannerChechikLangbergVertexSimple : public SpannerModule<TWeight>
{
public:
    SpannerChechikLangbergVertexSimple()
        : m_r(0)
        , m_seed(-1)
    {
    }

    void setFaultSize(int r)
    {
        m_r = r;
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
        double integralPart;
        if (std::modf(stretch, &integralPart) != 0.0)
        {
            error = "The stretch is required to be an odd integer, not " + to_string(stretch);
            return false;
        }
        int intStretch = static_cast<int>(stretch);
        if (intStretch < 1)
        {
            error = "The stretch must be >= 1.0";
            return false;
        }
        if (std::modf((stretch + 1) / 2, &integralPart) != 0.0)
        {
            error = "The stretch is required to be an odd integer, not " + to_string(stretch);
            return false;
        }
        if (GA.has(GraphAttributes::edgeDoubleWeight) || GA.has(GraphAttributes::edgeIntWeight))
        {
            error = "The graph must be unweighted";
            return false;
        }
        return true;
    }

private:
    const Graph *m_G;  //!< const reference to the original graph

    // the parameter k derived from the stretch
    int m_k;
    // the fault-tolerance parameter
    int m_r;
    // seed for random operations
    int m_seed;
    // calculated probability from m_r and number of nodes
    double m_prob;

    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_k = (stretch + 1) / 2;
        m_G = &GA.constGraph();

        int n = m_G->numberOfNodes();
        m_prob = pow(m_r + 3.0, 1.0 / (m_k - 1.0)) * pow(n / log(n), -1.0 / m_k);

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
        assertTimeLeft();

        // initialize a[i]
        std::vector<std::vector<node>> a(m_k + 1);
        a[0].resize(m_G->numberOfNodes());
        int count = 0;
        for (const node v : m_G->nodes)
        {
            a[0][count++] = v;
        }

        for (int i = 1; i <= m_k - 1; i++)
        {
            for (const node n : a[i - 1])
            {
                double randomD = randomDouble(0.0, 1.0);
                if (randomD <= m_prob)
                {
                    a[i].push_back(n);
                }
            }
        }

        std::vector<NodeArray<TWeight>> delta;
        // init delta (distances of nodes to each of the witnesses)
        delta.resize(m_k + 1);

        for (int t = 1; t <= m_r; t++)
        {
            // get all subset of nodes F of size t
            std::vector<std::vector<node>> node_subset_t;
            std::vector<int> combination;
            Array<node> all_nodes_for_combination;
            m_G->allNodes(all_nodes_for_combination);
            combination.resize(t);
            getSubsetPossibilities(t, 0, all_nodes_for_combination, combination, node_subset_t);

            // for every subset of nodes F of size t
            for (auto F : node_subset_t)
            {
                // looks expensive to do inside this loop
                // create subgraph without node set F
                GraphCopySimple graph_without_F(*m_G);

                EdgeArray<TWeight> m_weights;
                m_weights.init(graph_without_F, 1);
                delta[m_k].init(graph_without_F, std::numeric_limits<TWeight>::max());

                // remove nodes from F of nodes for cluster / spanner calculation
                for (node n : F)
                {
                    graph_without_F.delNode(graph_without_F.copy(n));
                }

                // cluster / spanner calculation for nodes without F
                // do k SSSP calls
                for (int i = m_k - 1; i >= 0; i--)
                {
                    assertTimeLeft();

                    node source = graph_without_F.newNode();

                    for (const node n : a[i])
                    {
                        // skip possibly deleted nodes in graph_without_F
                        if (graph_without_F.copy(n) == nullptr)
                        {
                            continue;
                        }
                        edge e = graph_without_F.newEdge(source, graph_without_F.copy(n));
                        m_weights[e] = 0.0;
                    }

                    // compute SSSP for newly inserted source vertex
                    NodeArray<edge> preds(graph_without_F);
                    NodeArray<TWeight> distResult(graph_without_F);
                    NodeArray<std::vector<node>> clusters(graph_without_F);
                    NodeArray<NodeArray<TWeight>> distanceReturns(graph_without_F);

                    // call dijkstra to get shortest distances
                    Dijkstra<TWeight>().call(graph_without_F, m_weights, source, preds, distResult,
                                             false);

                    delta[i] = distResult;

                    graph_without_F.delNode(source);

                    for (const node w : a[i])
                    {
                        if (graph_without_F.copy(w) == nullptr)
                        {
                            continue;
                        }
                        if (std::find(a[i + 1].begin(), a[i + 1].end(), graph_without_F.copy(w)) ==
                            a[i + 1].end())
                        {
                            modDijkstra(graph_without_F, graph_without_F.copy(w),
                                        clusters[graph_without_F.copy(w)], delta,
                                        distanceReturns[graph_without_F.copy(w)], i);
                        }
                    }
                }
            }
        }

        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
	 * @brief calculate all possible node combinations of size t
	 *
	 * @param t size of a node combination
	 * @param offset helping variable to determine combinations
	 * @param nodes array of all nodes of a graph
	 * @param current_combination helping vector to store combinations across each function call
	 * @param node_subset_t std::vector<std::vector<<node>> to store node combinations in
	 */
    void getSubsetPossibilities(int t, int offset, Array<node> &nodes,
                                std::vector<int> &current_combination,
                                std::vector<std::vector<node>> &node_subset_t)
    {
        if (t == 0)
        {
            // at this point, another combination is finished
            return;
        }
        for (int i = offset; i <= m_G->numberOfNodes() - t; i++)
        {
            // only get combinations with increasing node indices to avoid duplicates
            int current_size = current_combination.size();
            if (current_size - t - 1 >= 0 and
                nodes[i]->index() <= current_combination[current_size - t - 1])
            {
                continue;
            }
            current_combination[current_size - t] = i;
            getSubsetPossibilities(t - 1, offset + 1, nodes, current_combination, node_subset_t);

            // if last number of new combination got set write combination in node_subset_t
            if (t == 1)
            {
                std::vector<node> new_combination;
                new_combination.resize(current_size);
                node_subset_t.push_back(new_combination);
                for (int f = 0; f < current_size; f++)
                {
                    node_subset_t[node_subset_t.size() - 1][f] = nodes[current_combination[f]];
                }
            }
        }
        return;
    }

    /**
	 * @brief Executes a modified dijkstra search to form clusters. The modification to the standard
	 * dijkstra algorithm is that we relax an edge (u,v) only if the new distance found for
	 * the vertex v is smaller then the \p delta value of v of the next \p iteration. In this
	 * application an additional constraint has been set to enforce the trimmed cluster definition.
	 *
	 * @param graph_without_F Subgraph of original graph without set of nodes F
	 * @param source The start node of the search
	 * @param cluster Contains all the nodes visited by the search
	 * @param delta The precalculated distances to the samples
	 * @param distanceReturn The distances from the \p source node
	 * @param iteration The current iteration
	 */
    void modDijkstra(const GraphCopySimple &graph_without_F, const node &source,
                     std::vector<node> &cluster, const std::vector<NodeArray<TWeight>> &delta,
                     NodeArray<TWeight> &distanceReturn, const int iteration)
    {
        assertTimeLeft();

        // run modified dijkstra
        PrioritizedMapQueue<node, TWeight, std::less<TWeight>> pq(graph_without_F);
        distanceReturn.init(graph_without_F, std::numeric_limits<TWeight>::max());

        for (const node v : graph_without_F.nodes)
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
            for (adjEntry adj : graph_without_F.original(u)->adjEntries)
            {
                edge e = graph_without_F.copy(adj->theEdge());
                if (e == nullptr)
                {
                    // original adjEntries may include edges removed in graph_without_F
                    continue;
                }
                node v = e->opposite(u);

                // check conditions for relaxation
                auto weight = 1;
                if (distanceReturn[v] > distanceReturn[u] + weight &&
                    delta[iteration + 1][v] > distanceReturn[u] + weight &&
                    m_k >= distanceReturn[u] + weight)
                {
                    distanceReturn[v] = distanceReturn[u] + weight;
                    pq.decrease(v, distanceReturn[v]);
                    if (m_spanner->searchEdge(m_spanner->copy(graph_without_F.original(u)),
                                              m_spanner->copy(graph_without_F.original(v))) ==
                        nullptr)
                    {
                        m_spanner->newEdge(graph_without_F.original(e));
                        (*m_inSpanner)[graph_without_F.original(e)] = true;
                    }
                }
            }
        }
    }

    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_stretch;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};

}  // namespace ogdf
