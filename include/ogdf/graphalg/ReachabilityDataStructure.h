/** \file
 * \brief Implementation of the data structure described by
 * G.F. Italiano in 1986.
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
#include <ogdf/basic/Array2D.h>
#include <ogdf/basic/simple_graph_alg.h>

namespace ogdf {

/**
 * Reachability data structure.
 * 
 * Italiano. Amortized Efficiency of a Path Retrieval Data Structure.
 * Dipartimento di Informatica e Sistemistica, Univerity of Rome. 1986
 * 
 * The data structure takes O(n^2) time to initialize and supports directed
 * edge insertions in O(n) time (amortized). It is possible to anwer a
 * reachability query in constant time (true or false question).
 */
class ReachabilityDataStructure
{
    struct STNode {
        node vertex;
        std::vector<int> children;
        STNode(){};
        STNode(node v)
            : vertex(v){};
    };
    std::vector<STNode> nodeSet;
    NodeArray<NodeArray<int>> index;
    const Graph *m_G;

public:
    /**
     * Initialize the data structure
     *
     * @param g The graph the structure refers to
     */
    void init(Graph &g)
    {
        m_G = &g;
        index.init(*m_G);
        nodeSet.resize(m_G->numberOfNodes());
        int i = 0;
        for (node n : m_G->nodes)
        {
            nodeSet[i].vertex = n;
            index[n].init(*m_G, -1);
            index[n][n] = i++;
        }
    }

    /**
     * Add the directed edge (i,j) into structure
     *
     * @param i Start node
     * @param j End node
     */
    void add(node i, node j)
    {
        if (index[i][j] == -1)
        {
            for (node x : m_G->nodes)
            {
                if (index[x][i] != -1 && index[x][j] == -1)
                {
                    meld(x, j, i, j);
                }
            }
        }
    }

    /**
     * Check if node \p j is reachable from \p i
     *
     * @param i Start node
     * @param j End node
     * @return True if j is reachable from i
     */
    bool reachable(node i, node j)
    {
        return index[i][j] != -1;
    }

private:
    /**
     * Merges tree of \p x and the subtree of \p j rooted at node \p v.
     * Note that the tree of j is not changed.
     * 
     * @param x Node of which the tree should be merged
     * @param j The descendants of this node could be reached
     * @param u The node in the tree of x to which the pruned subtree will be linked
     * @param v Gives root of subtree of j
     */
    void meld(node x, node j, node u, node v)
    {
        nodeSet.emplace_back(v);
        index[x][v] = nodeSet.size() - 1;
        nodeSet[index[x][u]].children.push_back(nodeSet.size() - 1);
        for (int childIndex : nodeSet[index[j][v]].children)
        {
            node childVertex = nodeSet[childIndex].vertex;
            if (index[x][childVertex] == -1)
            {
                meld(x, j, v, childVertex);
            }
        }
    }
};
}  // namespace ogdf
