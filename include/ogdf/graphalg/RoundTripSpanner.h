/** \file
 * \brief Base class for all roundtrip spanner algorithms
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
 * Base class for roundtrip spanner algorithms. Contains static method to 
 * check for correctness and some functions used by multiple roundtrip
 * spanner algorithms.
 * 
 * @ingroup ga-spanner
 * @tparam TWeight
 */
template <typename TWeight>
class RoundTripSpanner : public SpannerModule<TWeight>
{
public:
    /**
	 * Validates a roundtrip spanner.
	 *
	 * @returns true iff the \p spanner is a roundtrip spanner of \p GA
	 */
    static bool isRoundTripSpanner(const GraphAttributes &GA, const GraphCopySimple &spanner,
                                   double stretch)
    {
        const Graph &origGraph = GA.constGraph();
        EpsilonTest eps;
        EdgeArray<double> origWeights(origGraph);
        EdgeArray<double> spannerWeights(spanner);
        NodeArray<edge> predecessors;
        NodeArray<NodeArray<double>> distanceMatrixOrig(origGraph);
        NodeArray<NodeArray<double>> distanceMatrixSpanner(spanner);
        for (edge e : origGraph.edges)
        {
            origWeights[e] = GA.doubleWeight(e);
        }
        for (edge e : spanner.edges)
        {
            spannerWeights[e] = GA.doubleWeight(spanner.original(e));
        }
        for (node n : origGraph.nodes)
        {
            node spannerN = spanner.copy(n);
            Dijkstra<double>().callUnbound(origGraph, origWeights, n, predecessors,
                                           distanceMatrixOrig[n], true, false);
            Dijkstra<double>().callUnbound(spanner, spannerWeights, spannerN, predecessors,
                                           distanceMatrixSpanner[spannerN], true, false);
        }
        for (node n : origGraph.nodes)
        {
            for (node m : origGraph.nodes)
            {
                if (n != m)
                {
                    double rtDistOrig = distanceMatrixOrig[n][m] + distanceMatrixOrig[m][n];
                    double rtDistSpanner = distanceMatrixSpanner[spanner.copy(n)][spanner.copy(m)] +
                                           distanceMatrixSpanner[spanner.copy(m)][spanner.copy(n)];
                    if (eps.greater(rtDistSpanner, stretch * rtDistOrig))
                    {
                        return false;
                    }
                }
            }
        }
        return true;
    }

protected:
    /**
	* Methdod returns the logarithm of parameter a to base b.
	*/
    double logbase(double a, double b)
    {
        return log(a) / log(b);
    }

    /**
     * Calculate an in tree or out tree and add edges of the tree to spanner.
     *
	 * @param graph The current graph
     * @param n The start node
     * @param preds The predecessor relations
     */
    void treeCalculation(const GraphCopySimple &graph, node n, const NodeArray<edge> &preds)
    {
        node curTarget = n;
        edge curEdge = preds[n];
        while (curEdge != nullptr)
        {
            edge origEdge = graph.original(curEdge);
            if (m_spanner->searchEdge(m_spanner->copy(origEdge->source()),
                                      m_spanner->copy(origEdge->target()), true) == nullptr)
            {
                edge newEdge = m_spanner->newEdge(origEdge);
                (*m_inSpanner)[origEdge] = true;
            }
            node pathNode = curEdge->opposite(curTarget);
            curEdge = preds[pathNode];
            curTarget = pathNode;
        }
    }

    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};
}  // namespace ogdf
