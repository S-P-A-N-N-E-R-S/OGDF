/** \file
 * \brief Implements the possibility to remove all degree two nodes and replace the
 * incident edges with one new edge. The weights of the two incident
 * edges are added to receive the weight of the new edge. If a shorter detour was found
 * this (smaller) value is set.
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

#include <ogdf/graphalg/GraphSimplification.h>

namespace ogdf {

GraphSimplification::GraphSimplification(const GraphAttributes &GA)
    : m_pGraph()
    , m_vOrig()
    , m_GA()
    , m_vSimplification()
    , m_eSimplification()
{
    m_pGraph = &GA.constGraph();

    Graph::construct(*m_pGraph, m_vSimplification, m_eSimplification);
    m_vOrig.init(*this);
    m_GA.init(*this, GA.attributes());

    for (node v : m_pGraph->nodes)
    {
        m_vOrig[m_vSimplification[v]] = v;
    }
    for (edge e : m_pGraph->edges)
    {
        if (m_GA.has(GraphAttributes::edgeDoubleWeight))
        {
            m_GA.doubleWeight(m_eSimplification[e]) = GA.doubleWeight(e);
        }
        else if (m_GA.has(GraphAttributes::edgeIntWeight))
        {
            m_GA.intWeight(m_eSimplification[e]) = GA.intWeight(e);
        }
    }
    NodeArray<bool> processed(*m_pGraph, false);
    Queue<node> queue;
    for (node n : this->nodes)
    {
        if (n->degree() == 2)
        {
            processed[m_vOrig[n]] = true;
            queue.append(n);
        }
    }
    while (!queue.empty())
    {
        vertexProcessing(queue, processed);
    }
}

void GraphSimplification::vertexProcessing(Queue<node> &queue, NodeArray<bool> &processed)
{
    node v = queue.pop();
    if (v->degree() == 2)
    {
        edge firstEdge = v->firstAdj()->theEdge();
        edge secondEdge = v->lastAdj()->theEdge();
        node neighbor1 = firstEdge->opposite(v);
        node neighbor2 = secondEdge->opposite(v);
        edge fEdge = this->searchEdge(neighbor1, neighbor2);
        this->delNode(v);
        if ((fEdge == nullptr) ||
            (m_GA.has(GraphAttributes::edgeDoubleWeight) &&
             m_GA.doubleWeight(fEdge) >
                 m_GA.doubleWeight(firstEdge) + m_GA.doubleWeight(secondEdge)) ||
            (m_GA.has(GraphAttributes::edgeIntWeight) &&
             m_GA.intWeight(fEdge) > m_GA.intWeight(firstEdge) + m_GA.intWeight(secondEdge)))
        {
            edge newEdge = this->newEdge(neighbor1, neighbor2);
            if (fEdge != nullptr)
            {
                this->delEdge(fEdge);
            }
            if (m_GA.has(GraphAttributes::edgeDoubleWeight))
            {
                m_GA.doubleWeight(newEdge) =
                    m_GA.doubleWeight(firstEdge) + m_GA.doubleWeight(secondEdge);
            }
            else if (m_GA.has(GraphAttributes::edgeIntWeight))
            {
                m_GA.intWeight(newEdge) = m_GA.intWeight(firstEdge) + m_GA.intWeight(secondEdge);
            }
        }
        if (neighbor1->degree() == 2 && !processed[m_vOrig[neighbor1]])
        {
            processed[m_vOrig[neighbor1]] = true;
            queue.append(neighbor1);
        }
        if (neighbor2->degree() == 2 && !processed[m_vOrig[neighbor2]])
        {
            processed[m_vOrig[neighbor2]] = true;
            queue.append(neighbor2);
        }
    }
}

}  // namespace ogdf