/** \file
 * \brief Implementation of the (2k+e) roundtrip spanner 
 * algorithm described by Roditty, Thorup and Zwick.
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
 * Multiplicative (2k+e) roundtrip spanner.
 *
 * Iam Roditty, Mikkel Thorup and Uri Zwick. Roundtrip Spanners and Roudtrip Routing in Directed Graphs.
 * ACM Transactions on Algorithms Vol 4, No 3, 2008
 *
 * Chun Jiang Zhu and Kam-Yiu Lam. Deterministic improved round-trip spanners. 
 * Information Processing Letters, 129:27-60. 2018
 *
 * Conditions for the graph:
 * - directed
 * - weighted
 *
 * The stretch value must be integer.
 * Technically the weights must be >=1. This is not checked though.
 *
 * The preconditions can be checked with RoundTripSpannerRodittyThorup::preconditionsOk.
 *
 * The (2k+e)-roundtrip-spanner has O((k^2/e)n^(1+1/k)log(nW)) edges, with
 * W beeing the maximum edge weight. In the paper no runtime complexity is proven.
 * 
 * The algorithm by Roditty was improved in 2018 by Zhu and Lam. They provided a derandomization of the cover
 * creation, which also leads to improved spanner sizes and runtime in practice. The size bound is improved by a factor of k.
 * These improvements are also implemented. The version can be set with the according method.
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class RoundTripSpannerRodittyThorup : public RoundTripSpanner<TWeight>
{
    // original graph
    const Graph *m_G;
    // the parameter k derived from the stretch
    int m_k;
    double m_epsilon;
    bool m_deterministic;
    // max edge weight
    TWeight m_W;
    EdgeArray<TWeight> m_weights;

public:
    RoundTripSpannerRodittyThorup(double epsilon = 1, bool deterministic = true)
        : RoundTripSpanner<TWeight>()
    {
        setEpsilon(epsilon);
        setDeterministic(deterministic);
    }

    void setDeterministic(bool setTo)
    {
        m_deterministic = setTo;
    }

    void setEpsilon(double epsilon)
    {
        OGDF_ASSERT(epsilon > 0);
        m_epsilon = epsilon;
    }

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
        if (intStretch <= 1)
        {
            error = "The stretch must be > 1.0";
            return false;
        }
        double integralPart;
        if (std::modf(stretch, &integralPart) != 0.0)
        {
            error = "The value of stretch has to be an integer";
            return false;
        }
        if (not(GA.has(GraphAttributes::edgeDoubleWeight) ||
                GA.has(GraphAttributes::edgeIntWeight)))
        {
            error = "The graph must be weighted";
            return false;
        }
        if (int((stretch - m_epsilon) / 2) == 0)
        {
            error = "Stretch and epsilon combination leads to k=0";
            return false;
        }
        if (m_epsilon <= 0)
        {
            error = "Epsilon must be > 0";
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
    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G = &GA.constGraph();
        m_k = (m_stretch - m_epsilon) / 2;
        // set weights and max edge weight
        m_weights.init(*m_G);
        m_W = 0;
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
        double epsilonStroke = m_epsilon / (2 * m_k);
        if (!m_deterministic)
        {
            for (int i = 1; i <= this->logbase(2 * m_G->numberOfNodes() * m_W, 1 + epsilonStroke);
                 i++)
            {
                assertTimeLeft();
                Array<std::vector<node>> rtCover;
                cover(pow(1 + epsilonStroke, i), rtCover);
            }
        }
        else
        {
            for (int i = 1; i <= this->logbase(2 * m_G->numberOfNodes() * m_W, 1 + epsilonStroke);
                 i++)
            {
                assertTimeLeft();
                std::vector<std::vector<node>> rtCover;
                coverDeterministic(pow(1 + epsilonStroke, i), rtCover);
            }
        }
        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    /**
	 * Creates a set of vertices, which are in a radius of maximum \p r (considering roundtrip distances). 
	 * The center of the this \p ball is node \p v.
	 *
	 * @param vertexSetV The vertices of the induced subgraph
	 * @param v The center of the ball
	 * @param r The radius
	 * @param ball All the nodes of the ball
	 * @param ballVertices Quick access to the ballVertices
	 */
    void ballCalculation(List<node> &vertexSetV, node v, TWeight r, std::vector<node> &ball,
                         NodeArray<bool> *ballVertices = nullptr)
    {
        // create subgraph of original graph induced by the set of vertices in vertexSetV
        GraphCopySimple inducedGraph;
        inducedSubGraph(*m_G, vertexSetV.begin(), inducedGraph);

        NodeArray<edge> predsIn;
        NodeArray<edge> predsOut;
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

        // only nodes of the vertexSetV can be in the ball
        for (node n : vertexSetV)
        {
            // check if roundtrip distance is smaller then radius
            if (distResult[inducedGraph.copy(n)] + distResult2[inducedGraph.copy(n)] <= r)
            {
                // insert node into ball
                ball.push_back(n);
                // only true for core ball
                if (ballVertices != nullptr)
                {
                    (*ballVertices)[n] = true;
                }
                // calculate shortest out paths and add edges to spanner
                this->treeCalculation(inducedGraph, inducedGraph.copy(n), predsIn);
                // calculate shortest in paths and add edges to spanner
                this->treeCalculation(inducedGraph, inducedGraph.copy(n), predsOut);
            }
        }
    }

    /**
	 * Creates a set of vertices, which are in a radius of maximum \p r (considering roundtrip distances). 
	 * The center of the this ball is node \p v.
	 *
	 * @param vertexSetV The vertices of the induced subgraph
	 * @param v The center of the ball
	 * @param r The radius
	 * @param ball All the nodes of the ball
	 * @param coreBallUnion The core ball
	 * @param ballVertices Quick access to the ball vertices
	 */
    void ballCalculationDeterministic(List<node> &vertexSetV, node v, TWeight r,
                                      std::vector<node> &ball, std::vector<node> &coreBallUnion,
                                      NodeArray<bool> &ballVertices)
    {
        int i = 1;
        ball.push_back(v);
        while (ball.size() > pow(m_G->numberOfNodes(), 1.f / m_k) * coreBallUnion.size())
        {
            for (node n : ball)
            {
                if (ballVertices[n] == false)
                {
                    ballVertices[n] = true;
                    coreBallUnion.push_back(n);
                }
            }
            ball.clear();
            ballCalculation(vertexSetV, v, i * r, ball);
            i++;
        }
    }

    /**
	 * Creates a roundtrip cover of the graph. The array \p roundTripCover contains
	 * balls. Each ball is represented by a vector. The graph is always a reference to the original graph.
	 * Derandomized version.
	 *
	 * @param r The radius 
	 * @param roundTripCover The roundtrip covers
	 */
    void coverDeterministic(TWeight r, std::vector<std::vector<node>> &roundTripCover)
    {
        Array<List<node>> vertexSet;
        int currSetIndex = 0;
        vertexSet.grow(1);
        // initialize set 0 with all graph nodes
        for (node n : m_G->nodes)
        {
            vertexSet[0].pushBack(n);
        }

        while (vertexSet[currSetIndex].size() > 0)
        {
            ListPure<node> pureList = vertexSet[currSetIndex].getListPure();
            // take first element to have a deterministic algo
            node u = pureList.front();
            std::vector<node> ball;
            std::vector<node> coreBallUnion;
            NodeArray<bool> coreBallVertices(*m_G, false);
            ballCalculationDeterministic(vertexSet[currSetIndex], u, r, ball, coreBallUnion,
                                         coreBallVertices);
            roundTripCover.push_back(ball);
            vertexSet.grow(1);
            // remove all vertices, which are in the core
            for (node v : vertexSet[currSetIndex])
            {
                if (!coreBallVertices[v])
                {
                    vertexSet[currSetIndex + 1].pushBack(v);
                }
            }
            currSetIndex++;
        }
    }

    /**
	 * Creates a roundtrip cover of the graph. The array \p roundTripCover contains
	 * balls. Each ball is represented by a vector. The graph is always a reference to the original graph.
	 * Randomized version.
	 *
	 * @param r The radius 
	 * @param roundTripCover The roundtrip covers
	 */
    void cover(TWeight r, Array<std::vector<node>> &roundTripCover)
    {
        Array<List<node>> vertexSet(m_k);

        // init k-1th vertexSet with all graph nodes
        for (node n : m_G->nodes)
        {
            vertexSet[m_k - 1].pushBack(n);
        }

        for (int i = m_k - 1; i >= 0; i--)
        {
            // sample vector s
            std::vector<node> s;
            double exp = (float)-i / m_k;
            sample(s, vertexSet[i], pow(m_G->numberOfNodes(), exp));
            std::vector<node> coreBallUnion;
            // we will calculate one ball for every node in the sample
            int rtCoverSize = 0;
            roundTripCover.grow(s.size());
            NodeArray<bool> coreBallVertices(*m_G, false);
            for (node v : s)
            {
                std::vector<node> ball;
                ballCalculation(vertexSet[i], v, (i + 1) * r, ball);
                roundTripCover[rtCoverSize++] = ball;
                // coreBall is not reset in the loop so all vertices are appended
                ballCalculation(vertexSet[i], v, i * r, coreBallUnion, &coreBallVertices);
            }
            // update the vertexSet for the next iteration
            for (node v : vertexSet[i])
            {
                if (!coreBallVertices[v])
                {
                    vertexSet[i - 1].pushBack(v);
                }
            }
        }
    }

    void sample(std::vector<node> &s, const List<node> &vertexSet, double prob)
    {
        for (node n : vertexSet)
        {
            if (prob > randomDouble(0, 1))
            {
                s.push_back(n);
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
