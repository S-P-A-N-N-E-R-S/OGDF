/** \file
 * \brief Implementation of the O(kloglog(n)) multiplicative 
 * roundtrip spanner algorithm described by Chechik, Liu and Rotem.
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

#include <ogdf/graphalg/RoundTripSpanner.h>

namespace ogdf {

/**
 * Multiplicative roundtrip-spanner.
 *
 * Shiri Chechik et al. Constant girth approximation for directed graphs in subquadratic time.
 * Proceedings of the 52nd Annual ACM SIGACT Symposium on Theory of Computing. 2020
 * 
 * Conditions for the graph:
 * - directed
 * - weighted
 *
 * THE STRETCH VALUE IS MISAPPROPRIATED TO BE USED AS THE PARAMETER K
 * This value has to be an integer > 0 
 *
 * The preconditions can be checked with RoundTripSpannerChechikLiu::preconditionsOk.
 *
 * The algorithm computes a O(k*loglog(n)) multiplicative roundtrip spanner with O(n^(1+1/k)) edges
 * with an additional logarithmic factor. The time is O(m^(1+1/k)) again with an additional
 * logarthmic factor. The approach is based on a recursive calculation of the roundtrip covers by
 * growing inballs and outballs around a center v at the same rate.
 *
 * Note that the runtimes depend on the maximum edge weight because the removal is not implemented.
 *
 * Another thing to notice is that the authors mention to take the union of the balls contained
 * in the roundtrip covers. Since this would lead to all edges beeing added, the in and out trees
 * of the contained balls are included into the spanner instead.
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class RoundTripSpannerChechikLiu : public RoundTripSpanner<TWeight>
{
    // original graph
    const Graph *m_G;
    int m_k;
    EdgeArray<TWeight> m_weights;
    double m_maxWeight;

public:
    //! @copydoc ogdf::SpannerModule::preconditionsOk
    virtual bool preconditionsOk(const GraphAttributes &GA, double stretch,
                                 std::string &error) override
    {
        const Graph &G = GA.constGraph();
        if (not GA.directed())
        {
            error = "The graph must be directed";
            return false;
        }
        double integralPart;
        if (std::modf(stretch, &integralPart) != 0.0)
        {
            error = "The value of stretch has to be an integer";
            return false;
        }
        int intStretch = static_cast<int>(stretch);
        if (intStretch < 1)
        {
            error = "The stretch must be >= 1.0";
            return false;
        }
        if (not(GA.has(GraphAttributes::edgeDoubleWeight) ||
                GA.has(GraphAttributes::edgeIntWeight)))
        {
            error = "The graph must be weighted";
            return false;
        }
        return true;
    }

private:
    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G = &GA.constGraph();
        // use the stretch value as the parameter k
        m_k = m_stretch;
        m_weights.init(*m_G);
        m_maxWeight = 0;
        for (edge e : m_G->edges)
        {
            m_weights[e] = SpannerModule<TWeight>::getWeight(GA, e);
            if (m_weights[e] > m_maxWeight)
            {
                m_maxWeight = m_weights[e];
            }
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        int origNumberOfNodes = m_G->numberOfNodes();
        for (int i = 0; i <= ceil(log2(origNumberOfNodes * m_maxWeight)); i++)
        {
            assertTimeLeft();
            NodeArray<bool> nodesInGBool(*m_G, true);
            cover(nodesInGBool, pow(2, i));
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
	 * Creates a O(kloglog(n), R) roundtrip cover and returns it. This means the balls of the cover
	 * have a \p radius of at must O(kloglog(n))*R.
	 * First the method creates an induced subgraph of the original graph. The nodes for this are given by the \p buildGraphNodes array. 
	 * The array references the original graph and contains booleans, defining the nodes, which should be
	 * contained in the induced subgraph. The roundtrip cover is a vector of balls and the balls are represended as graph copies of
	 * the original graph.
	 *
	 * @param buildGraphNodes The nodes used to build the incuded subgraph
	 * @param radius The radius of the cover
	 * @return The roundtrip cover
	 */
    std::vector<GraphCopySimple> cover(NodeArray<bool> &buildGraphNodes, TWeight radius)
    {
        if (radius == 0)
        {
            std::vector<GraphCopySimple> toReturn;
            return toReturn;
        }
        // this graph is GraphCopySimple of m_G
        GraphCopySimple graph;
        // check if there whould at least be one node in the subgraph
        bool addedOne = false;
        List<node> nodeList;
        for (node n : m_G->nodes)
        {
            if (buildGraphNodes[n])
            {
                addedOne = true;
                nodeList.pushBack(n);
            }
        }
        if (addedOne == false)
        {
            std::vector<GraphCopySimple> toReturn;
            return toReturn;
        }
        // create induced subgraph given by buildGraphNodes array (which is converted into the nodeList)
        // Important: If called with a GraphCopy (like here), one can access mapping with
        // original graph (.copy/.original). This possibility is disabled if called with explicit mapping arrays
        inducedSubGraph(*m_G, nodeList.begin(), graph);

        int iIn = 0;
        int iOut = 0;
        TWeight r = 5 * m_k * radius * log(log(graph.numberOfNodes()));
        // if you don't set this to 0 if < 0. The result may not be a RT-Spanner
        if (r <= 0)
        {
            r = 0;
        }
        // seems to work better with random nodes
        node v = graph.chooseNode();

        // return statement will trigger eventually
        // the anchor is set above: At some point the buildGraphNodes array will only
        // contain false values, so the induced subgraph will be empty. This triggers an empty return.
        while (true)
        {
            // calculate in and out balls with (i+1)R
            GraphCopySimple inBall;
            GraphCopySimple outBall;
            NodeArray<edge> preds;
            NodeArray<TWeight> distResult;
            NodeArray<bool> insertedNodeMarksIn(graph, false);
            ballCalculation(graph, (iIn + 1) * radius, inBall, v, insertedNodeMarksIn, preds,
                            distResult, true);
            NodeArray<bool> insertedNodeMarksOut(graph, false);
            ballCalculation(graph, (iOut + 1) * radius, outBall, v, insertedNodeMarksOut, preds,
                            distResult, false);
            // check size condition
            if (std::min(inBall.numberOfNodes(), outBall.numberOfNodes()) >=
                float(3 * graph.numberOfNodes()) / 4)
            {
                // calculate ball with 2r+R radius
                GraphCopySimple ball;
                ballComplete(graph, 2 * r + radius, ball, v);
                //remove all vertices, which are both in outBall and inBall
                for (node n : graph.nodes)
                {
                    if (insertedNodeMarksIn[n] && insertedNodeMarksOut[n])
                    {
                        // access original because graph is a subgraph of m_G
                        buildGraphNodes[graph.original(n)] = false;
                    }
                }
                std::vector<GraphCopySimple> res = cover(buildGraphNodes, radius);
                res.push_back(ball);
                return res;
            }
            // calculate in ball with i*R
            GraphCopySimple inBall2;
            insertedNodeMarksIn.init(graph, false);
            ballCalculation(graph, iIn * radius, inBall2, v, insertedNodeMarksIn, preds, distResult,
                            true);
            // call goodCut
            if (goodCut(graph, inBall2, inBall))
            {
                // reinit because it should only contain nodes of inBall
                NodeArray<bool> newNodeSet(*m_G, false);
                // set inBall nodes to true
                for (node n : inBall.nodes)
                {
                    newNodeSet[inBall.original(n)] = true;
                }
                std::vector<GraphCopySimple> rtCover1 = cover(newNodeSet, radius);
                // set inBall2 nodes to false
                for (node n : inBall2.nodes)
                {
                    buildGraphNodes[inBall2.original(n)] = false;
                }
                std::vector<GraphCopySimple> rtCover2 = cover(buildGraphNodes, radius);
                // union
                for (const GraphCopySimple &gcBall : rtCover2)
                {
                    rtCover1.push_back(gcBall);
                }
                return rtCover1;
            }

            // calculate out ball with i*R
            GraphCopySimple outBall2;
            insertedNodeMarksOut.init(graph, false);
            ballCalculation(graph, iOut * radius, outBall2, v, insertedNodeMarksOut, preds,
                            distResult, false);
            // call goodCut
            if (goodCut(graph, outBall2, outBall))
            {
                // reinit because it should only contain nodes of inBall
                NodeArray<bool> newNodeSet(*m_G, false);
                // set outBall nodes to true
                for (node n : outBall.nodes)
                {
                    newNodeSet[outBall.original(n)] = true;
                }
                std::vector<GraphCopySimple> rtCover1 = cover(newNodeSet, radius);
                // set outBall2 nodes to false
                for (node n : outBall2.nodes)
                {
                    buildGraphNodes[outBall2.original(n)] = false;
                }
                std::vector<GraphCopySimple> rtCover2 = cover(buildGraphNodes, radius);
                // union
                for (const GraphCopySimple &gcBall : rtCover2)
                {
                    rtCover1.push_back(gcBall);
                }
                return rtCover1;
            }
            if (inBall2.numberOfEdges() <= outBall2.numberOfEdges() ||
                outBall2.numberOfNodes() >= float(3 * graph.numberOfNodes()) / 4)
            {
                iIn++;
            }
            else
            {
                iOut++;
            }
        }
    }

    /**
	 * Determines whether recursion on \p ball2 and then deleting the \p ball1 nodes is good progress
	 *
	 * @param G The current graph
	 * @param ball1 The first ball
	 * @param ball2 The second ball
	 * @return Whether the check in the method returned true
	 */
    bool goodCut(const Graph &G, const GraphCopySimple &ball1, const GraphCopySimple &ball2)
    {
        if (ball2.numberOfNodes() <= (3.f / 4) * G.numberOfNodes() &&
            ball2.numberOfNodes() <= pow(ball1.numberOfNodes(), float(m_k - 1) / m_k) *
                                         pow(G.numberOfNodes(), 1.f / m_k) &&
            ball2.numberOfEdges() <=
                std::max(((1 + 1.f / m_k) * ball1.numberOfEdges()),
                         float(pow(ball1.numberOfEdges(), float(m_k - 1) / m_k) *
                               pow(G.numberOfEdges(), 1.f / m_k))))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
	 * Creates a complete \p ball, which is the combination of the in ball and the out ball with a specified \p radius
	 * The provided graph is a copy of the original graph.
	 *
	 * @param graph The current graph
	 * @param raidus The radius of the complete ball
	 * @param ball The complete ball
	 * @param v The center of the in and out balls
	 */
    void ballComplete(const GraphCopySimple &graph, TWeight radius, GraphCopySimple &ball, node v)
    {
        List<node> ballVertices;
        if (radius == 0)
        {
            ballVertices.pushBack(graph.original(v));
            inducedSubGraph(*m_G, ballVertices.begin(), ball);
            return;
        }
        NodeArray<bool> insertedNodeMarksIn(graph, false);
        NodeArray<bool> insertedNodeMarksOut(graph, false);
        NodeArray<edge> preds1;
        NodeArray<edge> preds2;
        NodeArray<TWeight> distResult1;
        NodeArray<TWeight> distResult2;
        // ball is reset at the end of the method
        ballCalculation(graph, radius, ball, v, insertedNodeMarksIn, preds1, distResult1, true);
        ballCalculation(graph, radius, ball, v, insertedNodeMarksOut, preds2, distResult2, false);
        for (node n : graph.nodes)
        {
            if (distResult1[n] + distResult2[n] <= radius)
            {
                ballVertices.pushBack(graph.original(n));
                // in
                this->treeCalculation(graph, n, preds1);
                // out
                this->treeCalculation(graph, n, preds2);
            }
        }
        // build induced subgraph
        inducedSubGraph(*m_G, ballVertices.begin(), ball);
    }

    /**
	 * Creates a set of vertices which are in a radius of maximum \p radius (considering roundtrip distances). 
     *
     * @param graph The current graph
     * @param radius The radius of the ball
     * @param ball The induced graph containing all the ball nodes \p v
	 * @param v The center of the ball
	 * @param insertedNodeMarks Quick access to all the nodes in the ball
	 * @param preds The predecessors for the shortest path
	 * @param distResult The distances calculated by the shortest path search
	 * @param reversed Whether to calculate an inball or outball
	 */
    void ballCalculation(const GraphCopySimple &graph, TWeight radius, GraphCopySimple &ball,
                         node v, NodeArray<bool> &insertedNodeMarks, NodeArray<edge> &preds,
                         NodeArray<TWeight> &distResult, bool reversed)
    {
        assertTimeLeft();
        List<node> ballVertices;
        if (radius == 0)
        {
            ballVertices.pushBack(graph.original(v));
            distResult.init(graph, std::numeric_limits<TWeight>::max());
            // only v is contained
            inducedSubGraph(*m_G, ballVertices.begin(), ball);
            return;
        }
        EdgeArray<TWeight> inducedWeights(graph);
        for (edge e : graph.edges)
        {
            inducedWeights[e] = m_weights[graph.original(e)];
        }
        Dijkstra<TWeight>().callBound(graph, inducedWeights, v, preds, distResult, true, reversed,
                                      nullptr, radius);
        for (node n : graph.nodes)
        {
            if (distResult[n] <= radius)
            {
                insertedNodeMarks[n] = true;
                ballVertices.pushBack(graph.original(n));
            }
        }
        inducedSubGraph(*m_G, ballVertices.begin(), ball);
    }
    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_stretch;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};
}  // namespace ogdf
