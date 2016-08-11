#ifndef CLASP_GRAPH_STATS_CALCULATOR_H
#define CLASP_GRAPH_STATS_CALCULATOR_H

#include <unordered_map>
#include <htd/TreeDecomposition.hpp>
#include <htd/IterativeImprovementTreeDecompositionAlgorithm.hpp>
#include <htd/TreeDecompositionAlgorithmFactory.hpp>
#include <exst/incidence_graph_stats.h>
#include <exst/dependency_graph_stats.h>
#include <exst/ExstTypes.h>

namespace exst
{
    /**
     * class used for graph stats
     */
    class GraphStatsCalculator
    {
    public:
        /**
         * Constructor of the GraphStatsCalculator class.
         * @param atomIds
         * @return
         */
        GraphStatsCalculator(std::unordered_map<uint32_t, uint32_t> &atomIds) :
                atomIds(atomIds), incidenceGraphStats(atomIds)
        {
        }

        ///calculator for statistics of dependency graph
        DependencyGraphStats dependencyGraphStats;

        ///mapping for atoms before and a
        std::unordered_map<uint32_t, uint32_t> &atomIds;

        ///calculator for statistics of incidence graph
        IncidenceGraphStats incidenceGraphStats;
    };
}

#endif //CLASP_GRAPH_STATS_CALCULATOR_H
