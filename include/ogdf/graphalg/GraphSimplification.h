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
#pragma once

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/Queue.h>

namespace ogdf {

class OGDF_EXPORT GraphSimplification : public Graph
{
protected:
    const Graph *m_pGraph;    // original graph
    NodeArray<node> m_vOrig;  // corresponding nodes in original graph
    GraphAttributes m_GA;
    NodeArray<node> m_vSimplification;
    EdgeArray<edge> m_eSimplification;

    GraphSimplification()
        : m_pGraph(nullptr)
    {
    }

public:
    // construction
    explicit GraphSimplification(const GraphAttributes &GA);
    virtual ~GraphSimplification()
    {
    }

    // returns original graph
    const Graph &original() const
    {
        return *m_pGraph;
    }
    // returns original node
    node original(node v) const
    {
        return m_vOrig[v];
    }
    // returns the GraphAttribues of the simplified graph
    GraphAttributes &getGraphAttributes()
    {
        return m_GA;
    }

private:
    void vertexProcessing(Queue<node> &queue, NodeArray<bool> &processed);
};
}  // namespace ogdf