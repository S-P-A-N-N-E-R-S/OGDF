/** \file
 * \brief Implementation of the (2k-1) multiplicative 
 * roundtrip spanner algorithm described by Cen, Duan and Gu
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
#include <queue>
#include <unordered_map>

namespace ogdf {

/**
 * Multiplicative (2k-1) roundtrip spanner.
 *
 * Ruoxu Cen et al. Roundtrip Spanners with (2k-1) Stretch.
 * arXiv preprint arXiv:1911.12411. 2020
 * 
 * Conditions for the graph:
 * - directed
 * - weighted
 *
 * The stretch value must be integer.
 *
 * The preconditions can be checked with RoundTripSpannerCenDuan::preconditionsOk.
 * 
 * The algorithm runs in to phases. First a preprocessing step is executed and then
 * the algorithm creates a set of roundtrip covers. The InOutTrees of balls are added to 
 * the spanner. The preprocessing adds O(n^1+1/k * log(n)) edges. 
 * The stretch is (2k-1) with O(kn^1+1/k * log(nW)) edges. The runtime is O(kmnlog(W)).
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class RoundTripSpannerCenDuan : public RoundTripSpanner<TWeight>
{
    // original graph
    const Graph *m_G;
    // the parameter k derived from the stretch
    int m_k;
    // maximum edge weight
    TWeight m_W;
    EdgeArray<TWeight> m_weights;

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
        if (intStretch % 2 == 0)
        {
            error = "The stretch must be odd";
            return false;
        }
        for (edge e : G.edges)
        {
            if (SpannerModule<TWeight>::getWeight(GA, e) < 1)
            {
                error = "Weights < 1 can lead to problems";
                return false;
            }
        }
        return true;
    }

private:
    struct NodeRadiusMatch {
        node n;
        TWeight radius;
    };
    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G = &GA.constGraph();
        m_k = (m_stretch + 1) / 2;
        m_W = 0;
        m_weights.init(*m_G);
        // read weights and max edge weight
        for (edge e : m_G->edges)
        {
            m_weights[e] = SpannerModule<TWeight>::getWeight(GA, e);
            if (m_W < m_weights[e])
            {
                m_W = m_weights[e];
            }
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        // this holds u and R(u) values
        std::vector<NodeRadiusMatch> R;
        preprocessing(R);
        // sort for cover method
        std::sort(R.begin(), R.end(), [](const NodeRadiusMatch &lhs, const NodeRadiusMatch &rhs) {
            return lhs.radius > rhs.radius;
        });
        double epsilon = 1.d / (2 * m_k - 2);

        for (int p = 0;
             p <= floor(this->logbase((2 * m_G->numberOfNodes() * m_W), 1 + epsilon)) + 1; p++)
        {
            List<node> allNodes;
            m_G->allNodes(allNodes);
            assertTimeLeft();
            cover(allNodes, pow(1 + epsilon, p), R);
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
	 * Creates a roundtrip cover of the graph. Since the cover is not used by the construction
     * any further the cover is not returned. The vector \p R is sorted and the \p nodeSet list
     * constains nodes of the original graph.
     *
     * @param nodeSet The set of nodes the cover is calculated on
     * @param L The step size
     * @param R Holds maximum length values for nodes
	 */
    void cover(List<node> &nodeSet, double L, const std::vector<NodeRadiusMatch> &R)
    {
        while (nodeSet.size() > 0)
        {
            assertTimeLeft();
            // find node with max length R, which is contained in nodeSet
            NodeRadiusMatch uMatch = R.front();
            int currPos = 0;
            while (!nodeSet.search(uMatch.n).valid())
            {
                uMatch = R[currPos++];
            }
            TWeight step = std::min(uMatch.radius / (m_k - 1), L);
            // calculate h
            int h = 1;
            while (true)
            {
                std::vector<node> ball;
                // not sure if it would be more efficient to do one
                // more ballCalculation and only add those edges with final h
                ballCalculation(nodeSet, uMatch.n, h * step, ball, true, true);
                if (ball.size() < pow(m_G->numberOfNodes(), (float)h / m_k))
                {
                    break;
                }
                else
                {
                    h++;
                };
            }
            // calculate not strictly smaller ball and remove the vertices from U
            std::vector<node> toDelete;
            ballCalculation(nodeSet, uMatch.n, (h - 1) * step, toDelete);
            // remove all the vertices in the ball from the nodeSet
            for (node n : toDelete)
            {
                auto findPos = nodeSet.search(n);
                if (findPos != nodeSet.end())
                {
                    nodeSet.del(findPos);
                }
            }
        }
    }

    /**
     * Preprocess the graph by creating a hitting set and adding
     * all the inwards and outwards directed shortest path trees to the spanner centered at vertices
     * inside the hitting set. 
     *
     * @param R Holds maximum length values for nodes
     */
    void preprocessing(std::vector<NodeRadiusMatch> &R)
    {
        Array<std::vector<node>> ballCollection(m_G->numberOfNodes());
        int counter = 0;
        for (node u : m_G->nodes)
        {
            assertTimeLeft();
            // sort distances by their roundtrip distance to u (increasing) to find R(u)
            std::vector<NodeRadiusMatch> toSortDistances;
            NodeArray<edge> preds;
            NodeArray<TWeight> distResult;
            NodeArray<TWeight> distResult2;
            Dijkstra<TWeight>().call(*m_G, m_weights, u, preds, distResult, true, false);
            Dijkstra<TWeight>().call(*m_G, m_weights, u, preds, distResult2, true, true);
            for (node n : m_G->nodes)
            {
                if (n != u)
                {
                    double rtDist = distResult[n] + distResult2[n];
                    if (rtDist == std::numeric_limits<TWeight>::infinity())
                    {
                        rtDist = std::numeric_limits<TWeight>::max();
                    }
                    toSortDistances.push_back(NodeRadiusMatch{n, rtDist});
                }
            }
            std::sort(toSortDistances.begin(), toSortDistances.end(),
                      [](const NodeRadiusMatch &lhs, const NodeRadiusMatch &rhs) {
                          return lhs.radius < rhs.radius;
                      });
            TWeight maxLengthR =
                toSortDistances[ceil(pow(m_G->numberOfNodes(), 1 - (1.f / m_k)))].radius;
            R.push_back(NodeRadiusMatch{u, maxLengthR});

            std::vector<node> ball;
            List<node> allNodes;
            m_G->allNodes(allNodes);
            ballCalculation(allNodes, u, maxLengthR, ball);
            ballCollection[counter++] = ball;
        }
        // get hitting set and calculate shortest path trees
        std::vector<node> hitSet;
        hitSet = hittingSet(ballCollection);
        for (node t : hitSet)
        {
            // calculate sp tree for node t (adds edges to spanner)
            shortestPathTree(t);
        }
    }

    /**
     * Hitting set algorithm. Returns the nodes in the hitting set.
     *
     * @param sets The node sets
     */
    std::vector<node> hittingSet(const Array<std::vector<node>> &sets)
    {
        std::vector<node> hittingSet;
        std::unordered_map<node, std::vector<int>> nodeSetMapping;
        NodeArray<int> currCounterMem(*m_G);
        // ERROR IN OGDF: Use of max ordering and decrease operation leads to wrong top element returns
        // Workaround for decrease operation with multimap implemented
        std::multimap<int, node, std::greater<int>> counterStruct;
        for (int i = 0; i < sets.size(); i++)
        {
            for (node n : sets[i])
            {
                nodeSetMapping[n].push_back(i);
            }
        }
        for (auto it = nodeSetMapping.begin(); it != nodeSetMapping.end(); it++)
        {
            node n = it->first;
            counterStruct.insert(std::pair<int, node>((it->second).size(), n));
            currCounterMem[it->first] = (it->second).size();
        }
        // initialization done
        Array<bool> covered(sets.size());
        covered.fill(false);
        while (true)
        {
            bool end = true;
            for (bool b : covered)
            {
                if (!b)
                {
                    end = false;
                }
            }
            if (end)
            {
                break;
            }
            node best = counterStruct.begin()->second;
            counterStruct.erase(counterStruct.begin());
            hittingSet.push_back(best);
            for (int i : nodeSetMapping[best])
            {
                if (!covered[i])
                {
                    covered[i] = true;
                    // decrement all node counts for set i
                    for (node contN : sets[i])
                    {
                        int currCountValue = currCounterMem[contN];
                        std::pair<std::multimap<int, node>::iterator,
                                  std::multimap<int, node>::iterator>
                            range = counterStruct.equal_range(currCountValue);
                        std::multimap<int, node>::iterator foundAt;
                        bool found = false;
                        // elements get deleted so check first
                        for (auto rangeElem = range.first; rangeElem != range.second; ++rangeElem)
                        {
                            if (rangeElem->second == contN)
                            {
                                foundAt = rangeElem;
                                found = true;
                            }
                        }
                        if (found)
                        {
                            counterStruct.erase(foundAt);
                            currCounterMem[contN] = currCounterMem[contN] - 1;
                            counterStruct.insert(std::pair<int, node>(currCountValue - 1, contN));
                        }
                    }
                }
            }
        }
        return hittingSet;
    }

    /**
     * Create inwards and outwards directed spt and add edges to spanner.
     *
     * @param t The start node
     */
    void shortestPathTree(node t)
    {
        NodeArray<edge> preds;
        NodeArray<TWeight> distResult;
        // outwards
        Dijkstra<TWeight>().call(*m_G, m_weights, t, preds, distResult, true, false);
        insertSpannerEdgesInPreprocessing(preds);
        // inwards
        Dijkstra<TWeight>().call(*m_G, m_weights, t, preds, distResult, true, true);
        insertSpannerEdgesInPreprocessing(preds);
    }

    /**
     * Inserts edges in the spanner based on the predecessors of
     * a dijkstra shortest path call.
     *
     * @param preds The predecessor relations
     */
    void insertSpannerEdgesInPreprocessing(const NodeArray<edge> &preds)
    {
        for (node n : m_G->nodes)
        {
            node curTarget = n;
            edge curEdge = preds[n];
            while (curEdge != nullptr)
            {
                if (m_spanner->searchEdge(m_spanner->copy(curEdge->source()),
                                          m_spanner->copy(curEdge->target()), true) == nullptr)
                {
                    m_spanner->newEdge(curEdge);
                    (*m_inSpanner)[curEdge] = true;
                }
                node pathNode = curEdge->opposite(curTarget);
                curEdge = preds[pathNode];
                curTarget = pathNode;
            }
        }
    }

    /**
	 * Creates a set of vertices which are in a radius of maximum \p radius (considering roundtrip distances). 
     *
     * @param vertexSetV The vertex set the induced subgraph is calculated on
     * @param v The center of the ball
     * @param r The radius of the ball
     * @param ball Contains all the nodes in the ball around \p v
     * @param strictlySmaller Whether the radius threshold is strictly smaller
     * @param insertEdges Whether the call should insert edges into the spanner
	 */
    void ballCalculation(List<node> &vertexSetV, node v, TWeight r, std::vector<node> &ball,
                         bool strictlySmaller = false, bool insertEdges = false)
    {
        GraphCopySimple inducedGraph;
        NodeArray<node> orig2New;
        EdgeArray<edge> edgesOrig2New;
        NodeArray<edge> predsIn;
        NodeArray<edge> predsOut;
        inducedSubGraph(*m_G, vertexSetV.begin(), inducedGraph);

        // start a dijkstra from node v (stopping at r)
        NodeArray<TWeight> distResult;
        NodeArray<TWeight> distResult2;
        EdgeArray<TWeight> inducedWeights(inducedGraph);

        for (edge e : inducedGraph.edges)
        {
            inducedWeights[e] = m_weights[inducedGraph.original(e)];
        }
        // called with r as max dist because if one way dist is larger already, rt dist can't be shorter if edge weights are positive
        Dijkstra<TWeight>().callBound(inducedGraph, inducedWeights, inducedGraph.copy(v), predsOut,
                                      distResult, true, false, nullptr, r);

        // call with reveresed arcs
        Dijkstra<TWeight>().callBound(inducedGraph, inducedWeights, inducedGraph.copy(v), predsIn,
                                      distResult2, true, true, nullptr, r);

        for (node n : vertexSetV)
        {
            node nCopy = inducedGraph.copy(n);
            if (strictlySmaller)
            {
                if (distResult[nCopy] + distResult2[nCopy] < r)
                {
                    ball.push_back(n);
                }
            }
            else
            {
                if (distResult[nCopy] + distResult2[nCopy] <= r)
                {
                    ball.push_back(n);
                }
            }
            if (insertEdges)
            {
                // calculate in trees and add edges to spanner
                this->treeCalculation(inducedGraph, nCopy, predsIn);
                // calculate out trees and add edges to spanner
                this->treeCalculation(inducedGraph, nCopy, predsOut);
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
