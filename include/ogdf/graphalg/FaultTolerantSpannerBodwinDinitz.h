/** \file
 * \brief Implementation of the fault tolerant spanner described by
 * Bodwin, Dinitz, Robelle
 *
 * \author Tim Hartmann
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

#include <ogdf/graphalg/ReachabilityDataStructure.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Fault tolerant spanner.
 *
 * Bodwin, Dinitz, Robelle. Optimal Vertex Fault-Tolerant Spanners in Polynomial Time.
 * SODA pages 2924--2938. 2021
 *
 * Conditions for the graph:
 * - undirected
 * - weighted
 *
 * The stretch must be an integer.
 *
 * The preconditions can be checked with FaultTolerantSpannerBodwinDinitz::preconditionsOk.
 *
 * The algorithm calculates an f-Vertex fault tolerant (2k-1) spanner with high probability.
 *
 * The authors did not make clear how to choose the threshold value m_tau, so the threshold value was set to 0.1 based on some
 * experimental executions with C=12. With this value combination the algorithm seems to achieves high accuracies even with a small number of nodes.
 * If the parameter f is increased above 1 or the number of nodes exceed very small amounts. The accuracy should be very close to 100%.
 * In addition the value can also be set with the according method. A small value of m_tau leads to higher probabilities
 * of getting a fault tolerant spanner, but the number of edges increases.
 * For higher probabilities the constant c can be increased at the expense of a higher runtime.
 * If the input graph is known to have a large amount of nodes or only f values above 1 will be used, the parameters should be to
 * get spanners with lower edge counts. Theoretically tau >= 1/(16e) and c>=(64\log2)/(f\log(n))+256 leads to a success probability of 1/n.
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class FaultTolerantSpannerBodwinDinitz : public SpannerModule<TWeight>
{
    // original graph
    const Graph *m_G;
    EdgeArray<TWeight> m_weights;
    int m_f = 1;
    int m_k;
    double m_tau = 0.1;
    int m_c = 12;

protected:
    EpsilonTest m_eps;

public:
    //! @copydoc ogdf::SpannerModule::preconditionsOk
    virtual bool preconditionsOk(const GraphAttributes &GA, double stretch, std::string &error)
    {
        if (GA.directed())
        {
            error = "The graph must be undirected";
            return false;
        }
        double integralPart;
        if (std::modf(stretch, &integralPart) != 0.0)
        {
            error = "The stretch is required to be an integer, not " + to_string(m_stretch);
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
            error = "The stretch must be odd";
            return false;
        }
        if (!(GA.has(GraphAttributes::edgeDoubleWeight) || GA.has(GraphAttributes::edgeIntWeight)))
        {
            error = "The graph must be weighted";
            return false;
        }
        return true;
    }

    void setFaultSize(int f)
    {
        m_f = f;
    }
    void setTau(double tau)
    {
        m_tau = tau;
    }
    void setC(int c)
    {
        m_c = c;
    }

private:
    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G = &GA.constGraph();
        m_weights.init(*m_G);
        m_k = (stretch + 1) / 2;
        for (edge e : m_G->edges)
        {
            m_weights[e] = SpannerModule<TWeight>::getWeight(GA, e);
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        Array<std::vector<node>> V(m_c * pow(m_f, 3) * log(m_G->numberOfNodes()) + 1);
        Array<Graph> layeredGraphs(m_c * pow(m_f, 3) * log(m_G->numberOfNodes()) + 1);
        Array<ReachabilityDataStructure> reachStructures(
            m_c * pow(m_f, 3) * log(m_G->numberOfNodes()) + 1);
        // contains the 2k additional nodes for every vertex of m_G
        Array<NodeArray<Array<node>>> layeredGraphNodeMatching(
            m_c * pow(m_f, 3) * log(m_G->numberOfNodes()) + 1);
        EdgeArray<std::vector<int>> L(*m_G);
        preprocessing(V, layeredGraphs, reachStructures, layeredGraphNodeMatching, L);

        Array<edge> edges(m_GA->constGraph().numberOfEdges());
        int count = 0;
        for (edge e : m_GA->constGraph().edges)
        {
            edges[count++] = e;
        }
        // sort the edges by weight (nondecreasing)
        edges.quicksort(EdgeWeightComparator(*m_GA));
        int counter = 0;
        for (edge e : edges)
        {
            assertTimeLeft();
            node u = e->source();
            node v = e->target();
            int pathsFound = 0;
            for (int i : L[e])
            {
                if (reachStructures[i].reachable(layeredGraphNodeMatching[i][u][0],
                                                 layeredGraphNodeMatching[i][v][2 * m_k - 1]))
                {
                    pathsFound++;
                }
            }
            double p = 1 - (double(pathsFound) / L[e].size());
            if (m_eps.geq(p, m_tau))
            {
                // add edge to spanner
                if (m_spanner->searchEdge(m_spanner->copy(u), m_spanner->copy(v)) == nullptr)
                {
                    m_spanner->newEdge(e);
                    (*m_inSpanner)[e] = true;
                }
                for (int i : L[e])
                {
                    // update data structure
                    for (int j = 1; j < 2 * m_k; j++)
                    {
                        reachStructures[i].add(layeredGraphNodeMatching[i][u][j - 1],
                                               layeredGraphNodeMatching[i][v][j]);
                        reachStructures[i].add(layeredGraphNodeMatching[i][v][j - 1],
                                               layeredGraphNodeMatching[i][u][j]);
                    }
                }
            }
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
	 * Preprocessing phase. Calculation of the random vertex sets and initialization of
	 * the reachability data structure for the subgraphs.
	 *
	 * @param V Sampled vertex sets
	 * @param layeredGraphs The created layered graphs
	 * @param reachStructures Calculated data structures
	 * @param layeredGraphNodeMatching Matching of original nodes to nodes of the layered graphs
	 * @param L Index set
	 */
    void preprocessing(Array<std::vector<node>> &V, Array<Graph> &layeredGraphs,
                       Array<ReachabilityDataStructure> &reachStructures,
                       Array<NodeArray<Array<node>>> &layeredGraphNodeMatching,
                       EdgeArray<std::vector<int>> &L)
    {
        NodeArray<std::vector<int>> LVertices(*m_G);
        for (int i = 0; i < m_c * pow(m_f, 3) * log(m_G->numberOfNodes()); i++)
        {
            layeredGraphNodeMatching[i].init(*m_G);
            assertTimeLeft();
            for (node v : m_G->nodes)
            {
                if (randomDouble(0, 1) < 1.0 / (2 * m_f))
                {
                    V[i].push_back(v);
                    LVertices[v].push_back(i);
                }
            }
            // create layered graph (H_i has no edges)
            for (node u : V[i])
            {
                layeredGraphNodeMatching[i][u].init(2 * m_k);
                for (int j = 0; j < 2 * m_k; j++)
                {
                    layeredGraphNodeMatching[i][u][j] = layeredGraphs[i].newNode();
                    if (j > 0)
                    {
                        layeredGraphs[i].newEdge(layeredGraphNodeMatching[i][u][j - 1],
                                                 layeredGraphNodeMatching[i][u][j]);
                    }
                }
            }
            // initialize reachability data structure
            reachStructures[i].init(layeredGraphs[i]);
            for (edge e : layeredGraphs[i].edges)
            {
                reachStructures[i].add(e->source(), e->target());
            }
        }

        // create the index sets L
        for (edge e : m_G->edges)
        {
            assertTimeLeft();
            node s = e->source();
            node t = e->target();
            int i = 0;
            int j = 0;
            int sizeLs = LVertices[s].size();
            int sizeLt = LVertices[t].size();
            while (i < sizeLs && j < sizeLt)
            {
                if (LVertices[s][i] < LVertices[t][j])
                {
                    i++;
                }
                else if (LVertices[t][j] < LVertices[s][i])
                {
                    j++;
                }
                else
                {
                    L[e].push_back(LVertices[t][j]);
                    i++;
                    j++;
                }
            }
        }
    }

    struct EdgeWeightComparator {
        EdgeWeightComparator(const GraphAttributes &GA)
            : m_GA(GA){};
        EdgeWeightComparator() = delete;
        EdgeWeightComparator(const EdgeWeightComparator &) = delete;
        EdgeWeightComparator &operator=(const EdgeWeightComparator &) = delete;

        bool less(edge a, edge b) const
        {
            return m_eps.less(getWeight(m_GA, a), getWeight(m_GA, b));
        }

    private:
        const GraphAttributes &m_GA;
        const EpsilonTest m_eps;
    };
    using SpannerModule<TWeight>::getWeight;
    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_stretch;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};
}  // namespace ogdf
