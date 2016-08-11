#ifndef CLASP_DEPENDENCY_GRAPH_STATS_H
#define CLASP_DEPENDENCY_GRAPH_STATS_H

#include <unordered_map>
#include <list>
#include <iostream>
#include <exst/ExstTypes.h>
#include <exst/dependency_graph_stats.h>

namespace exst
{
    /**
     * class used for calculating and saving stats of the dependency graph and the actual dependency graph
     */
    class DependencyGraphStats
    {
    public:
        /**
         * adds a rule to the dependency graph
         * @param body body of the rule
         * @param head head of the rule
         */
        void addRuleDependencyGraph(std::list<lit_type> body, std::list<lit_type> head);

        /**
         * prints the number of edges and nodes in the dependency graph
         */
        void printDepGraph();

        /**
         * Getter for the Dependency Graph.
         * @return the dependency graph
         */
        MyGraph &getDependencyGraph()
        {
            return dependencyGraph;
        }

        uint32_t getEdgeCount()
        {
            return edgecount;
        }

    private:
        //dependency graph
        MyGraph dependencyGraph;
        //mapping from graph vertices to literals
        std::unordered_map<uint32_t, uint32_t> atomVertexMap;
        //
        MyGraph edgeMap;
        //number of edges in the dependency graph
        uint32_t edgecount = 0;
        //number of nodes in the dependency graph
        uint32_t nodecount = 0;
    };
}

#endif //CLASP_DEPENDENCY_GRAPH_STATS_H
