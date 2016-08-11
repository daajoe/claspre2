#ifndef CLASP_EXSTTYPES_H
#define CLASP_EXSTTYPES_H

#include <unordered_map>
#include <htd/LibraryInstance.hpp>
#include <htd/TreeDecomposition.hpp>
#include <htd/IterativeImprovementTreeDecompositionAlgorithm.hpp>
#include <htd/TreeDecompositionAlgorithmFactory.hpp>
#include <htd/IMutableMultiGraph.hpp>

namespace exst
{
    /**
     * possible types of edges in the exst graph
     */
    enum EdgeType
    {
        POS = 1, NEG = -1, HEAD = 2, BODY = 3
    };

    /**
     * enum to represent the sign of a literal
     */
    enum Sign
    {
        POSITIVE = 1, NEGATIVE = -1
    };

    /**
     * type to represent a literal in exst
     */
    struct lit_type
    {
        /**
         * Constructor for lit_type with id and sign.
         * @param id id of the literal
         * @param s sign of the literal
         * @return a new instance of lit_type
         */
        lit_type(uint32_t id, Sign s) : id(id), s(s)
        {
        }

        /**
         * Empty construtor for lit_type.
         * @return a new instance of lit_type
         */
        lit_type()
        {
        }

        /// id of the literal
        uint32_t id;
        /// sign of the literal
        Sign s;

        /**
         * overload == operator for lit_type
         */
        bool operator==(lit_type l)
        {
            return (id == l.id) && (s == l.s);
        }
    };

    /**
     * datatype for internal exst graph representation
     */
    typedef std::unordered_map<uint32_t, std::unordered_map<uint32_t, EdgeType>> MyGraph;

    /**
     * Returns the max value of the given values.
     * @param values a map of values to get the max value from
     * @return the max value in the map
     */
    uint32_t maxValue(std::unordered_map<uint32_t, uint32_t> values);

    /**
     * Returns the min value of the given values.
     * @param values a map of values to get the min value from
     * @return the min value in the map
     */
    uint32_t minValue(std::unordered_map<uint32_t, uint32_t> values);

    /**
     * Calculates the number of edges in the given graph.
     * @param graph the graph to get the edgecount from
     * @return the edgecount of the grup
     */
    uint32_t edgeCount(MyGraph &graph);

    /**
     * Generates a copy of the given graph.
     * @param graph the graph to copy
     * @return a copy of the Graph
     */
    MyGraph copyMyGraph(MyGraph graph);

    /**
     * Generates the string representation of the graph in DIMACS format.
     * @param graph the graph to get the DIMACS format from
     * @return the graph in DIMACS format
     */
    std::string getDIMACS(MyGraph &graph);

    /**
     * Generates the string representation of the graph in Gr format.
     * @param graph the graph to get the Gr format from
     * @return the graph in Gr Format
     */
    std::string getGrFormat(MyGraph &graph);
}

#endif //CLASP_EXSTTYPES_H