/** \file
 * \brief GraphML parser for the frontend .graphml files
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
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

namespace ogdf {

template <typename TWeight>
class FEGraphMLParser
{
public:
    /**
	 * Read file and write into graph and graph attributes
	 */
    static void readGraphML(GraphAttributes &GA, Graph &G, std::string path,
                            int costFunctionIndex = 0)
    {
        std::unordered_map<int, node> findNode;
        std::ifstream graphmlFile(path);
        std::string line;
        edge currEdge;
        std::string costFunctionStrID = "weight_" + std::to_string(costFunctionIndex);
        while (std::getline(graphmlFile, line))
        {
            if (line.find("<node id=") != std::string::npos)
            {
                std::size_t pos = line.find("=");
                std::string graphmlIDStr = std::string(line.begin() + pos + 2, line.end() - 2);
                int graphmlID = std::stoi(graphmlIDStr);
                findNode[graphmlID] = G.newNode();
            }
            else if (line.find("<edge id=") != std::string::npos)
            {
                std::size_t pos1 = line.find("source=");
                std::size_t pos2 = line.find("target=");
                std::string nodeID1Str =
                    std::string(line.begin() + pos1 + 8, line.begin() + pos2 - 2);
                std::string nodeID2Str = std::string(line.begin() + pos2 + 8, line.end() - 2);
                int node1ID = std::stoi(nodeID1Str);
                int node2ID = std::stoi(nodeID2Str);
                currEdge = G.newEdge(findNode[node1ID], findNode[node2ID]);
            }
            else if (line.find(costFunctionStrID) != std::string::npos)
            {
                std::size_t pos1 = line.find(costFunctionStrID) + 10;
                std::size_t pos2 = line.find("</data");
                if (pos2 != std::string::npos)
                {
                    std::string weightStr = std::string(line.begin() + pos1, line.begin() + pos2);
                    TWeight weight = std::stod(weightStr);
                    writeWeight(GA, currEdge, weight);
                }
            }
        }
    }

private:
    static void writeWeight(GraphAttributes &GA, edge e, TWeight weight);
};

template <>
inline void FEGraphMLParser<int>::writeWeight(GraphAttributes &GA, edge e, int weight)
{
    if (GA.has(GraphAttributes::edgeIntWeight))
    {
        GA.intWeight(e) = weight;
    }
    else
    {
        GA.intWeight(e) = 1;
    }
}

template <>
inline void FEGraphMLParser<double>::writeWeight(GraphAttributes &GA, edge e, double weight)
{
    if (GA.has(GraphAttributes::edgeDoubleWeight))
    {
        GA.doubleWeight(e) = weight;
    }
    else
    {
        GA.doubleWeight(e) = 1;
    }
}
}  // namespace ogdf
