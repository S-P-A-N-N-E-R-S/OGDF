/** \file
 * \brief Implementation of the resilient spanner described by
 * Ausiello, Franciosa and Italiano
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

#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Resilient spanner.
 *
 * Ausiello, Franciosa, Italiano. On Resilient Graph Spanners.
 * Algorithmica volume 74, 2015
 *
 * Conditions for the graph:
 * - undirected
 * - weighted
 * - 2-edge-connected
 *
 * The stretch must be an integer.
 *
 * The preconditions can be checked with ResilientSpanner::preconditionsOk.
 *
 * The algorithm creates a resilient spanner, which is tolerant to transient
 * edge failures. This is done by adding a set of edges to an existing
 * spanner. The resulting spanner has O(W*n^(3/2)) edges with W=W_max/W_min and runs
 * in time O(mn + n^2*log(n))
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class ResilientSpanner : public SpannerModule<TWeight>
{
    // original graph
    const Graph *m_G;
    EdgeArray<TWeight> m_weights;
    double m_sigma;

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
        if (!(GA.has(GraphAttributes::edgeDoubleWeight) || GA.has(GraphAttributes::edgeIntWeight)))
        {
            error = "The graph must be weighted";
            return false;
        }
        if (!isSimple(GA.constGraph()))
        {
            error = "Graph must be simple";
            return false;
        }
        if (!isTwoEdgeConnected(GA.constGraph()))
        {
            error = "Graph must be 2-edge-connected";
            return false;
        }
        return true;
    }

    /**
     * The the sigma value of the algorithm. This value
     * has to be greater or equal to the stretch.
     *
     * @param sigma The sigma value
     * @return Whether sigma was set or the default was chosen
     */
    bool setSigma(double sigma)
    {
        if (sigma >= m_stretch)
        {
            m_sigma = sigma;
            return true;
        }
        else
        {
            m_sigma = m_stretch;
            return false;
        }
    }

private:
    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        m_GA = &GA;
        m_stretch = stretch;
        // this is already a valid spanner
        m_spanner = &spanner;
        m_inSpanner = &inSpanner;
        m_G = &GA.constGraph();
        m_weights.init(*m_G);
        for (edge e : m_G->edges)
        {
            m_weights[e] = SpannerModule<TWeight>::getWeight(GA, e);
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        // simplified check that the input is already a valid spanner
        if (&(m_spanner->original()) != m_G || m_spanner->numberOfEdges() == 0)
        {
            return SpannerModule<TWeight>::ReturnType::Error;
        }
        EdgeArray<TWeight> inducedWeights(*m_spanner);
        for (node r : m_G->nodes)
        {
            assertTimeLeft();
            // it seems not to be possible to set the weight of a newly inserted edge directly
            // this is why all the edge weights have to be set again in each iteration
            for (edge e : m_spanner->edges)
            {
                inducedWeights[e] = m_weights[m_spanner->original(e)];
            }
            // fragile edges belong to spanner
            std::vector<edge> fragileEdges;
            node rInSpanner = m_spanner->copy(r);
            // find fragile edges
            List<edge> adjEdgesOfR;
            rInSpanner->adjEdges(adjEdgesOfR);
            for (edge e : adjEdgesOfR)
            {
                // calc fragility of e in spanner
                TWeight edgeWeight = inducedWeights[e];
                node oppNode = e->opposite(rInSpanner);
                NodeArray<edge> preds;
                NodeArray<TWeight> distResult1;
                NodeArray<TWeight> distResult2;
                // del edge e by setting the weight to inf
                inducedWeights[e] = std::numeric_limits<TWeight>::max();
                Dijkstra<TWeight>().callBound(*m_spanner, inducedWeights, rInSpanner, preds,
                                              distResult1, false, false, oppNode,
                                              std::numeric_limits<TWeight>::max());
                // reset the weight
                inducedWeights[e] = edgeWeight;

                Dijkstra<TWeight>().callBound(*m_spanner, inducedWeights, rInSpanner, preds,
                                              distResult2, false, false, oppNode,
                                              std::numeric_limits<TWeight>::max());

                double fragilitySpanner =
                    (double(distResult1[oppNode] - distResult2[oppNode]) / distResult2[oppNode]) +
                    1;

                // calc fragility in original graph
                // del edge e by setting the weight to inf
                m_weights[m_spanner->original(e)] = std::numeric_limits<TWeight>::max();
                Dijkstra<TWeight>().callBound(*m_G, m_weights, r, preds, distResult1, false, false,
                                              m_spanner->original(oppNode),
                                              std::numeric_limits<TWeight>::max());
                // reset weight (is the same as in the spanner)
                m_weights[m_spanner->original(e)] = edgeWeight;

                Dijkstra<TWeight>().callBound(*m_G, m_weights, r, preds, distResult2, false, false,
                                              m_spanner->original(oppNode),
                                              std::numeric_limits<TWeight>::max());

                double fragilityOrig = ((distResult1[m_spanner->original(oppNode)] -
                                         distResult2[m_spanner->original(oppNode)]) /
                                        distResult2[m_spanner->original(oppNode)]) +
                                       1;
                if (m_eps.greater(fragilitySpanner, std::max(m_sigma, fragilityOrig)))
                {
                    // e belongs to spanner
                    fragileEdges.push_back(e);
                }
            }
            // add parsimonious cycles
            if (fragileEdges.size() > 0)
            {
                parsimoniousCycles(r, fragileEdges);
            }
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
     * Method creates parsimonious cyles for every given fragile edge. These backup paths
     * are added directly into the spanner.
     *
     * @param r The root vertex
     * @param fragileEdges The fragile edges calculated before
     */
    void parsimoniousCycles(node r, const std::vector<edge> &fragileEdges)
    {
        // r belongs to the original graph
        NodeArray<edge> spt;
        NodeArray<node> apex(*m_G);
        // number of edges from G without S in the paths
        NodeArray<int> k(*m_G, 0);
        EdgeArray<bool> edgesInSPT(*m_G, false);
        // distances of paths
        NodeArray<TWeight> delta(*m_G, 0);
        // compute parimonious shortest path tree
        // spt belongs to the original graph
        parsimoniousSPT(r, spt);
        // compute apex, edgesInSPT, k and delta arrays
        for (node n : m_G->nodes)
        {
            if (n != r)
            {
                node curTarget = n;
                edge curEdge = spt[n];
                while (curEdge != nullptr)
                {
                    edgesInSPT[curEdge] = true;
                    delta[n] = delta[n] + m_weights[curEdge];
                    if (!(*m_inSpanner)[curEdge])
                    {
                        k[n]++;
                    }
                    node pathNode = curEdge->opposite(curTarget);
                    if (pathNode == r)
                    {
                        apex[n] = curTarget;
                    }
                    curEdge = spt[pathNode];
                    curTarget = pathNode;
                }
            }
        }
        // find best paths (smallest weight with smallest edge count)
        NodeArray<std::vector<edge>> bestPaths(*m_G);
        NodeArray<int> bestPathLength(*m_G, std::numeric_limits<int>::max());
        NodeArray<TWeight> bestPathWeight(*m_G, std::numeric_limits<TWeight>::max());
        for (node y : m_G->nodes)
        {
            if (y != r)
            {
                List<edge> adjEdgesOfY;
                y->adjEdges(adjEdgesOfY);
                for (edge e : adjEdgesOfY)
                {
                    // y has an adjacent edge which is not in the spt, so it can be used to change branches
                    if (!edgesInSPT[e])
                    {
                        node z = e->opposite(y);
                        // its possible that z == r
                        if (z != r)
                        {
                            // if the graph is 2-edge-connected the edge must be contained because of the apex definition
                            edge rApexZ = m_G->searchEdge(r, apex[z]);
                            if (apex[y] != apex[z] &&
                                std::find(fragileEdges.begin(), fragileEdges.end(),
                                          m_spanner->copy(rApexZ)) != fragileEdges.end())
                            {
                                std::vector<edge> gamma;
                                // insert sp between r and y
                                node curTarget = y;
                                edge curEdge = spt[y];
                                while (curEdge != nullptr)
                                {
                                    gamma.push_back(curEdge);
                                    node pathNode = curEdge->opposite(curTarget);
                                    curEdge = spt[pathNode];
                                    curTarget = pathNode;
                                }
                                // insert edge (y,z)
                                gamma.push_back(e);
                                // insert sp between z and apex(z)
                                curTarget = z;
                                curEdge = spt[z];
                                while (curEdge != nullptr)
                                {
                                    gamma.push_back(curEdge);
                                    node pathNode = curEdge->opposite(curTarget);
                                    curEdge = spt[pathNode];
                                    curTarget = pathNode;
                                }
                                // remove the edge from apex[z] to r (this edge is fragile)
                                gamma.pop_back();
                                int gammaPathLength = k[y] + k[z];
                                if (!(*m_inSpanner)[e])
                                {
                                    gammaPathLength++;
                                }
                                // in k[z] rApexZ can be counted or not
                                // it is counted if it is not contained in the spanner
                                if (!(*m_inSpanner)[rApexZ])
                                {
                                    gammaPathLength--;
                                }

                                TWeight weightOfGammaPath =
                                    delta[y] + m_weights[e] + (delta[z] - m_weights[rApexZ]);
                                if (m_eps.less(weightOfGammaPath, bestPathWeight[apex[z]]))
                                {
                                    bestPaths[apex[z]] = gamma;
                                    bestPathLength[apex[z]] = gammaPathLength;
                                    bestPathWeight[apex[z]] = weightOfGammaPath;
                                }
                                else if (m_eps.equal(weightOfGammaPath, bestPathWeight[apex[z]]) &&
                                         m_eps.less(gammaPathLength, bestPathLength[apex[z]]))
                                {
                                    bestPaths[apex[z]] = gamma;
                                    bestPathLength[apex[z]] = gammaPathLength;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (edge eSpanner : fragileEdges)
        {
            edge e = m_spanner->original(eSpanner);
            node v = e->opposite(r);
            // edge (r,v) is already in the spanner
            for (edge spEdge : bestPaths[v])
            {
                node y = spEdge->nodes()[0];
                node z = spEdge->nodes()[1];
                if (m_spanner->searchEdge(m_spanner->copy(y), m_spanner->copy(z)) == nullptr)
                {
                    edge newEdge = m_spanner->newEdge(spEdge);
                    (*m_inSpanner)[spEdge] = true;
                }
            }
        }
    }

    /**
     * Computes the parsimonous shortest path tree. This is the shortest path three in which
     * every shortest path has the smallest number of edges possible.
     *
     * @param r The root vertex
     * @param spt The calculated shortest path tree
     */
    void parsimoniousSPT(node &r, NodeArray<edge> &spt)
    {
        PrioritizedMapQueue<node, TWeight, std::less<TWeight>> pq(*m_G);
        NodeArray<TWeight> distances;
        distances.init(*m_G, std::numeric_limits<TWeight>::max());
        NodeArray<int> pathLength(*m_G, std::numeric_limits<int>::max());

        for (node v : m_G->nodes)
        {
            pq.push(v, distances[v]);
        }
        spt.init(*m_G, nullptr);
        pq.decrease(r, (distances[r] = 0));
        pathLength[r] = 0;
        while (!pq.empty())
        {
            node v = pq.topElement();
            pq.pop();
            if (distances[v] == std::numeric_limits<TWeight>::max())
            {
                break;
            }
            for (adjEntry adj : v->adjEntries)
            {
                edge e = adj->theEdge();
                node w = e->opposite(v);
                auto weight = m_weights[e];
                if (m_eps.greater(distances[w], distances[v] + weight))
                {
                    spt[w] = e;
                    pq.decrease(w, (distances[w] = distances[v] + weight));
                    if (!(*m_inSpanner)[e])
                    {
                        pathLength[w] = pathLength[v] + 1;
                    }
                }
                // Additional condition in comparision to standard dijkstra
                else if (m_eps.equal(distances[w], distances[v] + weight) &&
                         m_eps.greater(pathLength[w], pathLength[v] + 1))
                {
                    spt[w] = e;
                    if (!(*m_inSpanner)[e])
                    {
                        pathLength[w] = pathLength[v] + 1;
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
