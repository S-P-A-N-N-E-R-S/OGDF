/** \file
 * \brief Implementation of the O(klog(n)) multiplicative roundtrip
 * spanner algorithm described by Pachocki, Roditty and Sidford.
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
 * Pachocki, Roditty and Sidford. Approximating Cycles in Directed Graphs: Fast Algorithms for Girth and Rountrip Spanners
 * Proceedings of the 29nd Annual ACM-SIAM Symposium on Discrete Algorithms. 2018
 * 
 * Conditions for the graph:
 * - directed
 * - weighted
 *
 * THE STRETCH VALUE IS MISAPPROPRIATED TO BE USED AS THE PARAMETER K
 * This value has to be an integer > 0
 *
 * The preconditions can be checked with RoundTripSpannerPachockiRoditty::preconditionsOk.
 *
 * The algorithm creates a O(klog(n)) multiplicative roundtrip spanner with (n^(1+1/k)*log(n)*log(nW)) edges
 * with a high probability. The runtime is O(mn^(1/k)*log^4(n)*log(nW)). The spanner is created by calculating
 * a set of roundtrip covers with increasing radii. The constant C can be adapted to achieve higher
 * propability bounds.
 *
 * Note that the estimateBallSizes method differs from the method presented in the paper. The method is 
 * based on sampling a set of vertices. The set size is calculated by a function which depends on the
 * number of vertices. This function is not applicable to graphs with small node counts because
 * the amount of sampled vertices would be higher then the number of vertices. Instead of this function, a 
 * different function is used if the vertex count is lower then 2000.
 *
 * The authors give a method to remove the dependency on the maximum edge length. This method is
 * not implemented. Instead the simple modification presented by Chechik and Liu was chosen.
 *
 * Another thing to notice is that the authors mention to take the union of the balls contained
 * in the roundtrip covers. Since this would lead to all edges beeing added, the in and out trees
 * of the contained balls are included into the spanner instead.
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class RoundTripSpannerPachockiRoditty : public RoundTripSpanner<TWeight>
{
    // original graph
    const Graph *m_G;
    int m_k;
    EdgeArray<TWeight> m_weights;
    double m_maxWeight;

    const double C = 2.0;

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
        // is set to log2 in paper
        for (int i = 0; i <= ceil(log2(m_G->numberOfNodes() * m_maxWeight)); i++)
        {
            assertTimeLeft();
            fastCover(pow(2, i));
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
     * Constructs roundtrip covers with radius \p R. A roundtrip cover
     * is essentially the union of probabilistic covers.
     * Nothing is returned because a return is not necessary for the
     * roundtrip spanner construction.
     *
     * @param R The radius of the roundtrip cover
     */
    void fastCover(TWeight R)
    {
        double r = 6 * R * m_k * log(m_G->numberOfNodes());
        // constant is set globally
        double endAt =
            C * ceil(pow(m_G->numberOfNodes(), 1.f / m_k)) * ceil(log(m_G->numberOfNodes()));
        // initial call to cover with original graph
        NodeArray<bool> initialNodes(*m_G, true);
        for (int i = 1; i <= endAt; i++)
        {
            probabilisticCover(initialNodes, r);
        }
    }

    /**
     * Calculates one roundtrip cover probabilistically. The cover is
     * not returned since this is not required for the spanner calculation.
     *
     * @param buildGraphNodes The nodes used to build the incuded subgraph
     * @param r The radius of the cover
     */
    void probabilisticCover(NodeArray<bool> buildGraphNodes, TWeight r)
    {
        assertTimeLeft();
        // induced subgraph of m_G
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
            return;
        }
        inducedSubGraph(*m_G, nodeList.begin(), graph);

        // estimate the ball sizes for every node
        NodeArray<std::pair<double, double>> estimatedSizes;
        estimatedSizes = estimateBallSizes(graph, C * r, 1.f / 8);

        // fill the vertex sets by analysis of the estimated sizes
        std::vector<node> S_out;
        std::vector<node> S_in;
        NodeArray<bool> outHash(graph, false);
        NodeArray<bool> inHash(graph, false);
        std::vector<node> intersection;
        for (node v : graph.nodes)
        {
            if (estimatedSizes[v].first >= 3.f / 4)
            {
                S_out.push_back(v);
                outHash[v] = true;
            }
            if (estimatedSizes[v].second >= 3.f / 4)
            {
                S_in.push_back(v);
                inHash[v] = true;
            }
            if (estimatedSizes[v].first >= 3.f / 4 && estimatedSizes[v].second >= 3.f / 4)
            {
                intersection.push_back(v);
            }
        }

        // check if intersection between S_out and S_in is not empty
        if (intersection.size() > 0)
        {
            // choose random vertex from the intersection
            node u = intersection[randomNumber(0, intersection.size() - 1)];
            // calculate in and out balls
            std::vector<node> outBall;
            NodeArray<bool> nodesContainedOut(graph, false);
            NodeArray<edge> preds;
            NodeArray<TWeight> distResult;
            ballCalculation(graph, C * r, outBall, u, preds, distResult, nodesContainedOut, false);
            std::vector<node> inBall;
            NodeArray<bool> nodesContainedIn(graph, false);
            ballCalculation(graph, C * r, inBall, u, preds, distResult, nodesContainedIn, true);
            // check size of intersection
            int intersectionSize = 0;
            for (node n : graph.nodes)
            {
                if (nodesContainedOut[n] && nodesContainedIn[n])
                {
                    intersectionSize++;
                }
            }
            if (intersectionSize < 1.f / 4 * graph.numberOfNodes())
            {
                return;
            }
            // pick rB and use it for the complete ball calculation
            double rB = randomDouble(2 * C * r, 2 * (C + 1) * r);
            std::vector<node> compBall;
            ballComplete(graph, rB, compBall, u);
            // update boolean array so the induced subgraph can be build in the recursive call
            // remove all the vertices which are in the complete ball
            for (node n : compBall)
            {
                buildGraphNodes[n] = false;
            }
            // recurse
            probabilisticCover(buildGraphNodes, r);
            return;
        }
        std::vector<std::vector<node>> clusters;
        if (S_out.size() <= 1.f / 2 * graph.numberOfNodes())
        {
            // calculate out clusters
            std::vector<node> vertexSet;
            for (node n : graph.nodes)
            {
                if (!outHash[n])
                {
                    vertexSet.push_back(n);
                }
            }
            clusterCalculation(graph, vertexSet, r, false, clusters);
        }
        else
        {
            // calculate in clusters
            std::vector<node> vertexSet;
            for (node n : graph.nodes)
            {
                if (!inHash[n])
                {
                    vertexSet.push_back(n);
                }
            }
            clusterCalculation(graph, vertexSet, r, true, clusters);
        }
        // find max size and check return criterium
        int maxSize = 0;
        for (std::vector<node> cluster : clusters)
        {
            if (cluster.size() > maxSize)
            {
                maxSize = cluster.size();
            }
        }
        if (maxSize > 7.f / 8 * graph.numberOfNodes())
        {
            return;
        }
        // prepare induced subgraphs for all the clusters
        for (int i = 0; i < clusters.size(); i++)
        {
            NodeArray<bool> inducedClusterNodes(*m_G, false);
            for (node n : clusters[i])
            {
                inducedClusterNodes[graph.original(n)] = true;
            }
            // recurse
            probabilisticCover(inducedClusterNodes, r);
        }
    }

    /**
	 * Creates a complete ball, which is the combination of the in ball and the out ball with a specified \p radius
	 * The provided graph is a copy of the original graph.
     *
     * @param graph The current graph
     * @param radius The radius of the complete ball
     * @param compBall All the nodes contained in the complete ball
     * @param v The center of the ball
	 */
    void ballComplete(const GraphCopySimple &graph, TWeight radius, std::vector<node> &compBall,
                      node v)
    {
        NodeArray<bool> nodesContainedHash(graph, false);
        NodeArray<edge> preds1;
        NodeArray<TWeight> distResult1;
        NodeArray<edge> preds2;
        NodeArray<TWeight> distResult2;
        std::vector<node> ball;
        // ball gets filled by these calls
        ballCalculation(graph, radius, ball, v, preds1, distResult1, nodesContainedHash, true);
        ballCalculation(graph, radius, ball, v, preds2, distResult2, nodesContainedHash, false);

        for (node n : graph.nodes)
        {
            if (distResult1[n] + distResult2[n] <= radius)
            {
                compBall.push_back(graph.original(n));
                // in
                this->treeCalculation(graph, n, preds1);
                // out
                this->treeCalculation(graph, n, preds2);
            }
        }
    }

    /**
	 * Creates a set of vertices which are in a radius of maximum \p radius (considering roundtrip distances). 
     *
     * @param graph The current graph
     * @param radius The radius of the ball
     * @param ball The nodes contained in the ball around center \p v
	 * @param v The center of the ball
	 * @param preds The predecessors for the shortest path
	 * @param distResult The distances calculated by the shortest path search
     * @param nodesContainedHash Quick access to the contained nodes
	 * @param reversed Whether to calculate an inball or outball
	 */
    void ballCalculation(const GraphCopySimple &graph, TWeight radius, std::vector<node> &ball,
                         node v, NodeArray<edge> &preds, NodeArray<TWeight> &distResult,
                         NodeArray<bool> &nodesContainedHash, bool reversed)
    {
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
                nodesContainedHash[n] = true;
                ball.push_back(graph.original(n));
            }
        }
    }

    /**
     * This method calculates for every node the estimated size of a ball by sampling a set of vertices
     * and then calculating the roundtrip distance from every node in the graph to the sampled vertices.
     * Note the difference to the presented function in the paper.
     *
     * @param graph The current graph
     * @param r The radius
     * @param epsilon Value used in the sample function
     * @return The estimated in and out ball sizes for each node
     */
    NodeArray<std::pair<double, double>> estimateBallSizes(const GraphCopySimple &graph, TWeight r,
                                                           double epsilon)
    {
        NodeArray<std::pair<double, double>> estimatedSizes(graph);
        EdgeArray<TWeight> inducedWeights(graph);
        for (edge e : graph.edges)
        {
            inducedWeights[e] = m_weights[graph.original(e)];
        }
        // PRESENTED FUNCTION IN PAPER NOT APPLICABLE FOR SMALL NUMBER OF VERTICES
        // USAGE OF DIFFERENT FUNCTION FOR VERTEX COUNTS SMALLER THEN 2000
        // The function was choosen in such a way that for small graphs a large percentage is sampled.
        // A smooth transition to the presented function is acceived by using a linear function.
        int t;
        if (graph.numberOfNodes() > 2000)
        {
            t = ceil(20 * log(graph.numberOfNodes() / pow(epsilon, 2)));
        }
        else if (graph.numberOfNodes() <= 100)
        {
            t = 4 * pow(graph.numberOfNodes(), 0.35);
            if (t > graph.numberOfNodes())
            {
                t = graph.numberOfNodes();
            }
        }
        else
        {
            t = 0.113158 * graph.numberOfNodes() + 8.68421;
        }
        Array<node> samples(t);
        Array<node> allNodesContainer;
        graph.allNodes(allNodesContainer);
        for (int i = 0; i < t; i++)
        {
            int nodeIndexSample = randomNumber(0, graph.numberOfNodes() - 1);
            samples[i] = allNodesContainer[nodeIndexSample];
        }

        for (node s : samples)
        {
            NodeArray<edge> preds;
            NodeArray<TWeight> distResult;
            // call dijkstra with max length r
            // then you can just add the reached vertices directly
            // out
            Dijkstra<TWeight>().callBound(graph, inducedWeights, s, preds, distResult, true, false,
                                          nullptr, r);

            // in
            NodeArray<TWeight> distResult2;
            Dijkstra<TWeight>().callBound(graph, inducedWeights, s, preds, distResult2, true, true,
                                          nullptr, r);
            // compute fraction of vertices which have been reached
            for (node v : graph.nodes)
            {
                if (distResult[v] <= r)
                {
                    estimatedSizes[v].first = estimatedSizes[v].first + 1.f / samples.size();
                }
                if (distResult2[v] <= r)
                {
                    estimatedSizes[v].second = estimatedSizes[v].second + 1.f / samples.size();
                }
            }
        }
        return estimatedSizes;
    }

    /**
     * The method uses the exponential distribution to assign vertices in the graph to
     * clusters. The given \p vertexSet contains the cluster centers.
     *
     * @param graph The current graph
     * @param vertexSet The cluster centers
     * @param r The radius
     * @param inCluster Whether to calculate an in or out cluster
     * @param clusters The resulting clusters
     */
    void clusterCalculation(const GraphCopySimple &graph, const std::vector<node> &vertexSet,
                            TWeight r, bool inCluster, std::vector<std::vector<node>> &clusters)
    {
        // induced weights for dijkstra call
        EdgeArray<TWeight> inducedWeights(graph);
        for (edge e : graph.edges)
        {
            inducedWeights[e] = m_weights[graph.original(e)];
        }
        // pick for every vertex in vertexSet a value from the exponential distribution
        double beta = log(graph.numberOfNodes()) / r;
        std::default_random_engine generator;
        std::exponential_distribution<double> distribution(beta);
        NodeArray<double> x(graph);
        for (node v : vertexSet)
        {
            x[v] = distribution(generator);
        }

        // init cluster structures
        // not assigned vertices are placed at the end
        clusters.resize(vertexSet.size() + 1);

        NodeArray<TWeight> minValues(graph, std::numeric_limits<TWeight>::max());
        // -1 means that the node is not assigned to any cluster
        NodeArray<int> clusterCenters(graph, -1);
        int clusterCounter = 0;
        // loop through possible clusters
        for (node v : vertexSet)
        {
            NodeArray<edge> preds;
            NodeArray<TWeight> distResult;
            Dijkstra<TWeight>().callBound(graph, inducedWeights, v, preds, distResult, true,
                                          inCluster,  // set if in or out cluster
                                          nullptr, std::numeric_limits<TWeight>::max());

            // update the current centers and min values for every node in graph
            for (node n : graph.nodes)
            {
                if (-x[v] + distResult[n] < minValues[n] && -x[v] + distResult[n] <= 0)
                {
                    clusterCenters[n] = clusterCounter;
                    minValues[n] = -x[v] + distResult[n];
                }
            }
            clusterCounter++;
        }

        // assign clusters
        for (node n : graph.nodes)
        {
            if (clusterCenters[n] != -1)
            {
                clusters[clusterCenters[n]].push_back(n);
            }
            else
            {
                // n was not assigned to any cluster
                clusters[vertexSet.size()].push_back(n);
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
