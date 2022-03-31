/** \file
 * \brief Implementation of a k-spanner exact algorithm by
 * Ahmed et al.
 *
 * \author Dennis Benz
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

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/external/coin.h>
#include <ogdf/graphalg/SpannerModule.h>

#include <iomanip>

namespace ogdf {

/**
 * Minimum spanner by solving a binary LP
 *
 * R. Ahmed, K. Hamm, M. Jebelli, S. Kobourov, F. Sahneh und R. Spence.
 * Approximation algorithms and an integer program for multi-level graph spanners.
 * International Symposium on Experimental Algorithms. Springer. 2019, S. 541–562.
 *
 * Conditions for the graph:
 *  - simple
 *  - undirected
 *  - connected
 *
 * The stretch \f$s\f$ must satisfy: \f$s\geq1\in\mathbb{R}\f$.
 *
 * The preconditions can be checked with SpannerAhmed::preconditionsOk.
 *
 * Calculates a minimum \f$k\f$-spanner by solving a polynomial sized ILP.
 * The ILP has \f$\mathcal{O}(|E||K|)\f$ variables and constraints with
 * \f$|K|=\mathcal{O}(|V|^2)\f$.
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class SpannerAhmed : public SpannerModule<TWeight>
{
public:
    SpannerAhmed()
        : m_osi(nullptr)
    {
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
        if (!isConnected(GA.constGraph()))
        {
            error = "The graph is not connected";
            return false;
        }
        if (stretch < 1.0)
        {
            error = "The stretch must be >= 1.0";
            return false;
        }
        return true;
    }

private:
    const Graph *m_G;      //!< const reference to the original graph
    GraphCopySimple m_rG;  //!< reduced graph from original graph

    // Solving the LP
    std::unique_ptr<OsiSolverInterface> m_osi;  //!< the solver
    std::vector<CoinPackedVector>
        m_constraints;  //!< Holds all constraints so they can be freed at destruction.

    /**
	 * Resets the LP defining variables.
	 */
    void resetLP()
    {
        m_constraints.clear();
    }

    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        resetLP();
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);

        m_osi = std::unique_ptr<OsiSolverInterface>(CoinManager::createCorrectOsiSolverInterface());
        m_osi->messageHandler()->setLogLevel(0);  // 0=nothing .. 4=verbose

        m_G = &GA.constGraph();
        m_rG.init(*m_G);
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        // Calculate all-pairs shortest paths
        EdgeArray<TWeight> edgeCosts(*m_G);
        for (edge weightEdge : m_G->edges)
        {
            edgeCosts[weightEdge] = getWeight(*m_GA, weightEdge);
        }
        NodeArray<NodeArray<TWeight>> shortestPathMatrix(*m_G);
        dijkstra_SPAP(*m_G, shortestPathMatrix, edgeCosts);
        assertTimeLeft();

        // Reduce graph by removing edges
        reduceGraph(shortestPathMatrix);
        assertTimeLeft();

        if (m_rG.numberOfNodes() == 0 || m_rG.numberOfEdges() == 0)
        {
            return SpannerModule<TWeight>::ReturnType::Feasible;
        }

        // Create graph copy with directed edges
        EdgeArray<std::pair<edge, edge>> directedEdges(m_rG);
        GraphCopySimple directedCopy;
        directedCopy.init(m_rG);
        directedCopy.reverseAllEdges();
        for (edge originalEdge : m_rG.edges)
        {
            edge reversedEdge = directedCopy.copy(originalEdge);
            edge undirectedEdge = directedCopy.newEdge(originalEdge);
            directedEdges[originalEdge] = {undirectedEdge, reversedEdge};
        }
        assertTimeLeft();

        // Set up ILP
        m_osi->setObjSense(1);  // minimize

        // One column per edge (x_e)
        CoinPackedVector zero;
        EdgeArray<int> edgeIndices(m_rG);
        int numVar = 0;
        for (edge e : m_rG.edges)
        {
            int lb = getEdgeLowerBound(m_rG.original(e), shortestPathMatrix);
            m_osi->addCol(zero, lb, 1,
                          getWeight(*m_GA, m_rG.original(e)));  // vec, lower, upper, objective
            m_osi->setInteger(numVar);
            edgeIndices[e] = numVar++;
        }
        assertTimeLeft();

        // Arrays for edge (i,j) on u-v-path in Spanner (x_((i,j))^(uv))
        NodeArray<NodeArray<EdgeArray<int>>> pathEdgeIndices(directedCopy);
        for (node u : directedCopy.nodes)
        {
            pathEdgeIndices[u].init(directedCopy);
            for (node v : directedCopy.nodes)
            {
                pathEdgeIndices[u][v].init(directedCopy);
            }
        }

        // One column for edge (i,j) on u-v-path in Spanner (x_((i,j))^(uv))
        for (node u : directedCopy.nodes)
        {
            for (node v : directedCopy.nodes)
            {
                // One direction
                if (u->index() < v->index())
                {
                    for (edge ij : directedCopy.edges)
                    {
                        m_osi->addCol(zero, 0, 1, 0);  // vec, lower, upper, objective
                        pathEdgeIndices[u][v][ij] = numVar++;
                    }
                }
            }
        }
        assertTimeLeft();

        // Create Constraints
        int numConst = 0;
        for (node u : directedCopy.nodes)
        {
            node original_u = m_rG.original(directedCopy.original(u));
            for (node v : directedCopy.nodes)
            {
                node original_v = m_rG.original(directedCopy.original(v));
                if (u->index() < v->index())
                {
                    // Spanner constraint <= t * d_G(u,v)
                    m_constraints.emplace_back();
                    CoinPackedVector &spannerConstraint = m_constraints.back();
                    for (edge ij : directedCopy.edges)
                    {
                        edge originalEdge_ij = m_rG.original(directedCopy.original(ij));
                        spannerConstraint.insert(pathEdgeIndices[u][v][ij],
                                                 getWeight(*m_GA, originalEdge_ij));
                    }
                    // Sense: 'E' ==   'G' >=   'L' <=
                    m_osi->addRow(spannerConstraint, 'L',
                                  m_stretch * shortestPathMatrix[original_u][original_v],
                                  0);  // constraint, sense, rhs (will be set below), range
                    numConst++;

                    // Multicommodity constrains / simple path constraints
                    for (node i : directedCopy.nodes)
                    {
                        // Flow constraint
                        m_constraints.emplace_back();
                        CoinPackedVector &flowConstraint = m_constraints.back();
                        // Outgoing edges
                        List<edge> outEdges;
                        i->outEdges(outEdges);
                        for (edge outEdge : outEdges)
                        {
                            flowConstraint.insert(pathEdgeIndices[u][v][outEdge], 1);
                        }
                        // Ingoing edges
                        List<edge> inEdges;
                        i->inEdges(inEdges);
                        for (edge inEdge : inEdges)
                        {
                            flowConstraint.insert(pathEdgeIndices[u][v][inEdge], -1);
                        }
                        int rhs;
                        if (i->index() == u->index())
                        {
                            rhs = 1;
                        }
                        else if (i->index() == v->index())
                        {
                            rhs = -1;
                        }
                        else
                        {
                            rhs = 0;
                        }
                        m_osi->addRow(flowConstraint, 'E', rhs, 0);
                        numConst++;

                        // Single outgoing edge constraint
                        m_constraints.emplace_back();
                        CoinPackedVector &singleEdgeConstraint = m_constraints.back();
                        for (edge outEdge : outEdges)
                        {
                            singleEdgeConstraint.insert(pathEdgeIndices[u][v][outEdge], 1);
                        }
                        m_osi->addRow(singleEdgeConstraint, 'L', 1, 0);
                        numConst++;
                    }

                    // Edge constraint
                    for (edge e : m_rG.edges)
                    {
                        m_constraints.emplace_back();
                        CoinPackedVector &edgeConstraint = m_constraints.back();
                        edgeConstraint.insert(pathEdgeIndices[u][v][directedEdges[e].first], 1);
                        edgeConstraint.insert(pathEdgeIndices[u][v][directedEdges[e].second], 1);
                        edgeConstraint.insert(edgeIndices[e], -1);
                        m_osi->addRow(edgeConstraint, 'L', 0, 0);
                        numConst++;
                    }
                    assertTimeLeft();
                }
            }
        }

        // Solve ILP
        m_osi->branchAndBound();
        assertTimeLeft();

        if (m_osi->isProvenOptimal())
        {
            const double *solution = m_osi->getColSolution();
            for (edge e : m_rG.edges)
            {
                if (solution[edgeIndices[e]])
                {
                    (*m_inSpanner)[m_rG.original(e)] = true;
                    m_spanner->newEdge(m_rG.original(e));
                }
            }
            return SpannerModule<TWeight>::ReturnType::Optimal;
        }
        return SpannerModule<TWeight>::ReturnType::Error;
    }

    /**
     * Create a reduced graph by removing not possible edges of the original graph.
     * Such edge would never be included in the spanner, as less costly paths
     * are existing between its two nodes.
     *
     * @param shortestPathMatrix All-pairs-shortest-path matrix of original graph
     */
    void reduceGraph(const NodeArray<NodeArray<TWeight>> &shortestPathMatrix)
    {
        for (edge e : m_G->edges)
        {
            // Shortest path test
            if (shortestPathMatrix[e->source()][e->target()] < getWeight(*m_GA, e))
            {
                m_rG.delEdge(m_rG.copy(e));
                continue;
            }

            // Remove not possible edge
            bool possibleEdge = false;
            for (node u : m_G->nodes)
            {
                for (node v : m_G->nodes)
                {
                    if (u->index() < v->index())
                    {
                        TWeight d_ui = shortestPathMatrix[u][e->source()];
                        TWeight d_jv = shortestPathMatrix[e->target()][v];
                        TWeight d_uv = shortestPathMatrix[u][v];
                        if (d_ui + getWeight(*m_GA, e) + d_jv <= m_stretch * d_uv)
                        {
                            possibleEdge = true;
                        }
                        // Check reversed edge
                        d_ui = shortestPathMatrix[u][e->target()];
                        d_jv = shortestPathMatrix[e->source()][v];
                        if (d_ui + getWeight(*m_GA, e) + d_jv <= m_stretch * d_uv)
                        {
                            possibleEdge = true;
                        }
                    }
                }
            }
            if (!possibleEdge)
            {
                m_rG.delEdge(m_rG.copy(e));
            }
        }
    }

    /**
     * Determines the lower bound of the edge variable and sets it to one
     * if the edge is definitely in the spanner.
     *
     * @param e Edge corresponding to variable
     * @param shortestPathMatrix All-pairs-shortest-path matrix of original graph
     * @return Lower bound of edge variable
     */
    int getEdgeLowerBound(edge e, const NodeArray<NodeArray<TWeight>> &shortestPathMatrix)
    {
        GraphCopySimple copy(*m_G);
        copy.delEdge(copy.copy(e));

        EdgeArray<TWeight> edgeCosts(copy);
        for (edge weightEdge : copy.edges)
        {
            edgeCosts[weightEdge] = getWeight(*m_GA, copy.original(weightEdge));
        }
        NodeArray<NodeArray<TWeight>> reducedShortestPathMatrix(copy);
        dijkstra_SPAP(copy, reducedShortestPathMatrix, edgeCosts);

        for (node u : m_G->nodes)
        {
            for (node v : m_G->nodes)
            {
                if (u->index() < v->index())
                {
                    if (reducedShortestPathMatrix[copy.copy(u)][copy.copy(v)] >
                        m_stretch * shortestPathMatrix[u][v])
                    {
                        return 1;
                    }
                }
            }
        }
        return 0;
    }

    using SpannerModule<TWeight>::getWeight;
    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_stretch;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};

}  // namespace ogdf