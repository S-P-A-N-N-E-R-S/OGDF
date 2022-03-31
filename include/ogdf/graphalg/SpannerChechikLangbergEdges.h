/** \file
 * \brief Implementation of an r-edge-fault-tolerant (2k-1)-spanner algorithm
 * from Chechik and Langberg
 *
 * \author Julian Pascal Wittker
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

#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Algorithm for calculating fault-tolerant-spanners.
 *
 * Chechik, S., Langberg, M., Peleg, D., & Roditty, L. (2010).
 * Fault tolerant spanners for general graphs.
 * SIAM Journal on Computing, 39(7), 3403-3423.
 *
 * Conditions for the graph:
 * - simple
 * - undirected
 * - weighted
 *
 * The stretch mus satisfy \f$k\geq1\f$.
 *
 * The preconditions can be checked with SpannerChechikLangbergEdges::preconditionsOk.
 *
 * Calculates a r-edge-fault-tolerant (2k-1)-spanner with \f$r\geq1\f$
 *
 * @ingroup ga-spanner <- TODO: ?
 */
template <typename TWeight, typename TSpanner>
class SpannerChechikLangbergEdges : public SpannerModule<TWeight>
{
    static_assert(std::is_base_of<SpannerModule<TWeight>, TSpanner>::value,
                  "TSpanner must be a derived class of SpannerModule<TWeight> in "
                  "SpannerChechikLangbergEdges<TWeight, TSpanner>.");

public:
    SpannerChechikLangbergEdges()
        : m_r(0)
    {
    }

    void setFaultSize(int r)
    {
        m_r = r;
    }

    //! @copydoc ogdf::SpannerModule::preconditionsOk
    virtual bool preconditionsOk(const GraphAttributes &GA, double stretch,
                                 std::string &error) override
    {
        if (!isSimple(GA.constGraph()))
        {
            error = "The graph is not simple";
            return false;
        }
        double integralPart;
        if (std::modf(stretch, &integralPart) != 0.0)
        {
            error = "The stretch is required to be an integer, not " + to_string(m_stretch);
            return false;
        }
        if (stretch < 3.0)
        {
            error = "The stretch must be >= 3.0";
            return false;
        }
        if (std::fmod(stretch, 2) == 0)
        {
            error = "The stretch must be uneven";
            return false;
        }
        if (GA.directed())
        {
            error = "The graph must be undirected";
            return false;
        }
        if (!(GA.has(GraphAttributes::edgeDoubleWeight) || GA.has(GraphAttributes::edgeIntWeight)))
        {
            error = "The graph must be weighted";
            return false;
        }
        return true;
    }

private:
    const Graph *m_G;  //!< const pointer to the original graph

    // the fault-tolerance parameter
    int m_r;

    //! @copydoc ogdf::SpannerModule::init
    virtual void init(const GraphAttributes &GA, double stretch, GraphCopySimple &spanner,
                      EdgeArray<bool> &inSpanner) override
    {
        SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
        m_G = &GA.constGraph();
    }

    //! @copydoc ogdf::SpannerModule::execute
    virtual typename SpannerModule<TWeight>::ReturnType execute() override
    {
        // initialize empty edge set for spanner
        std::vector<edge> e_sp;
        GraphCopySimple copy_G(*m_G);
        GraphAttributes copy_ga(copy_G, m_GA->attributes());
        copy_ga.directed() = false;
        for (edge e : m_G->edges)
        {
            copy_ga.doubleWeight(copy_G.copy(e)) = m_GA->doubleWeight(e);
        }

        // for r iterations
        for (int i = 0; i < m_r + 1; i++)
        {
            assertTimeLeft();

            // remove edges from graph copy for next spanner algorithm
            for (auto edge : e_sp)
            {
                copy_ga.doubleWeight(edge) = std::numeric_limits<TWeight>::max();
            }
            // clear vector to only keep track of newly inserted edges
            e_sp.clear();

            GraphCopySimple spanner_copy_G(copy_G);
            EdgeArray<bool> e_Arr(copy_G);

            TSpanner copy_G_Spanner;
            if (this->isTimelimitEnabled())
            {
                copy_G_Spanner.setTimelimit(this->getTimeLeft());
            }
            auto res = copy_G_Spanner.call(copy_ga, m_stretch, spanner_copy_G, e_Arr);

            // adapt spanner's edge array and add edges to spanner's edge set
            for (edge e : copy_G.edges)
            {
                if (e_Arr[e])
                {
                    node a = e->nodes()[0];
                    node b = e->nodes()[1];
                    if (m_spanner->searchEdge(m_spanner->copy(copy_G.original(a)),
                                              m_spanner->copy(copy_G.original(b))) == nullptr)
                    {
                        m_spanner->newEdge(copy_G.original(e));
                        (*m_inSpanner)[copy_G.original(e)] = true;
                    }
                    e_sp.push_back(e);
                }
            }
        }

        return SpannerModule<TWeight>::ReturnType::Feasible;
    }

    using SpannerModule<TWeight>::getWeight;
    using SpannerModule<TWeight>::assertTimeLeft;
    using SpannerModule<TWeight>::m_GA;
    using SpannerModule<TWeight>::m_stretch;
    using SpannerModule<TWeight>::m_spanner;
    using SpannerModule<TWeight>::m_inSpanner;
};

}  // namespace ogdf
