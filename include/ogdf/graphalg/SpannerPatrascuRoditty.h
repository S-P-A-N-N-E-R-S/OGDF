/** \file
 * \brief Implementation of a (2d+1)-spanner algorithm
 * from Patrascu and Roditty
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
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/SpannerModule.h>
#include <algorithm>

#include <math.h>
#include <ogdf/basic/Array2D.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/graphalg/Dijkstra.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>

namespace ogdf {

/**
 * Algorithm for calculating (2d+1)-spanners.
 *
 * M. Patrascu and L. Roditty.
 * Distance Oracles beyond the Thorup-Zwick Bound (2010).
 * In Proceedings of the 2010 IEEE 51st Annual Symposium on Foundations of Computer Science (FOCS '10).
 * IEEE Computer Society, USA, 815–823.
 * DOI:https://doi.org/10.1109/FOCS.2010.83
 *
 * Conditions for the graph:
 * - simple
 * - undirected
 * - unweighted
 *
 * The stretch must be 2.
 *
 * The preconditions can be checked with SpannerPatrascuRoditty::preconditionsOk.
 *
 * Calculates a (2d+1)-spanner
 *
 * @ingroup ga-spanner
 */
template <typename TWeight>
class SpannerPatrascuRoditty : public SpannerModule<TWeight>
{
public:
    //! Initializes a spanner module.
    SpannerPatrascuRoditty()
        : m_seed(-1)
    {
    }

    //! @copydoc ogdf::SpannerModule::preconditionsOk
    bool preconditionsOk(const GraphAttributes &GA, double stretch, std::string &error) override
    {
        if (GA.directed())
        {
            error = "The graph must be undirected";
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
    // seed for random operations
    int m_seed;

    // original graph
    GraphCopySimple m_G;

    //! @copydoc ogdf::SpannerModule::init
    void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
              EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G.init(GA.constGraph());

        if (m_seed == -1)
        {
            setSeed(time(NULL));
        }
    }

    //! @copydoc ogdf::SpannerModule::execute
    typename SpannerModule<TWeight>::ReturnType execute() override
    {
        assertTimeLeft();

        double probabilityA = pow(m_G.numberOfNodes(), (-1.f / 3));
        double probabilityB = pow(m_G.numberOfNodes(), (-2.f / 3));

        std::vector<node> setA;
        std::vector<node> setB;
        std::vector<node> setC;

        NodeArray<std::vector<node>> balls(m_G);

        // create set A and B
        for (const node n : m_G.nodes)
        {
            double randomD = randomDouble(0, 1);
            if (randomD <= probabilityA)
            {
                setA.push_back(n);
            }
            randomD = randomDouble(0, 1);
            if (randomD <= probabilityB)
            {
                setB.push_back(n);
            }
        }

        // insert setA nodes in setC
        for (const node n : setA)
        {
            setC.push_back(n);
        }

        // grow a ball around each vertex of setB
        // and stop if a vertex of setA is reached
        for (const node n : setB)
        {
            growBall(n, setA, setC);
        }

        // grow a ball for the nodes in setC too.
        for (const node n : m_G.nodes)
        {
            growBall(n, setC, balls[n]);
        }

        return SpannerPatrascuRoditty::ReturnType::Feasible;
    }

    /**
	 * Grows a ball around node \p source and stop if a node
	 * of \p set is reached
	 *
	 * @param source The start node on which the ball starts to grow
	 * @param set Contains nodes which once reached by the ball stop the algorithm
	 * @param reachedNodes Set of nodes not in \p set which are reached by the ball around \p source
	 */
    void growBall(const node source, const std::vector<node> &set, std::vector<node> &reachedNodes)
    {
        assertTimeLeft();

        // run modified BFS
        PrioritizedMapQueue<node, int, std::less<int>> pq(m_G);

        NodeArray<int> distToReachedNodes(m_G, std::numeric_limits<int>::max());

        NodeArray<bool> visited(m_G, false);

        visited[source] = true;
        distToReachedNodes[source] = 0;
        pq.push(source, distToReachedNodes[source]);

        while (!pq.empty())
        {
            const node current = pq.topElement();
            pq.pop();

            if (std::find(set.begin(), set.end(), current) != set.end())
            {
                break;
            }

            if (std::find(reachedNodes.begin(), reachedNodes.end(), current) ==
                    reachedNodes.end() &&
                current != source)
            {
                reachedNodes.push_back(current);
            }

            for (const adjEntry adj : current->adjEntries)
            {
                const node w = adj->twinNode();
                if (visited[w] == false)
                {
                    visited[w] = true;
                    distToReachedNodes[w] = distToReachedNodes[current] + 1;
                    pq.push(w, distToReachedNodes[w]);

                    if (m_spanner->searchEdge(m_spanner->copy(m_G.original(current)),
                                              m_spanner->copy(m_G.original(w))) == nullptr)
                    {
                        edge e = adj->theEdge();
                        m_spanner->newEdge(m_G.original(e));
                        (*m_inSpanner)[m_G.original(e)] = true;
                    }
                }
            }
        }
    }

    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};
}  // namespace ogdf
