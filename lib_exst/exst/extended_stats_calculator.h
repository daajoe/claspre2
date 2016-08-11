#ifndef CLASP_EXTENDED_STATS_CALCULATOR_H
#define CLASP_EXTENDED_STATS_CALCULATOR_H

#include <unordered_map>
#include <list>
#include <exst/ExstTypes.h>
#include <exst/graph_stats_calculator.h>

namespace exst
{
    /**
     * class used to generate the extended stats
     */
    class StatsCalculator
    {
    public:

        /**
         * returns the single instance of the GraphStatsCalculator
         */
        static StatsCalculator &getInstance()
        {
            static StatsCalculator calc;
            return calc;
        }

        /**
         * adds a dependency to the dependency graph
         */
        void parseRule(std::list<lit_type> body, std::list<lit_type> head);

        /**
         * Maps the id of an atom before and after preprocessing.
         * @param before id of the atom before preprocessing
         * @param after id of the atom after preprocesing
         */
        void addId(uint32_t before, uint32_t after);

        /**
         * Used to set the Symbol Table.
         * @param table the symbol table of the program
         */
        void setSymbolTable(std::unordered_map<uint32_t, const char *> &table);

        /**
         * Prints the Statistics.
         */
        void printExtendedStats();

        /**
         * Starts the calculation of the stats, implicitly called by printExtendedStats.
         */
        void calculateStats();

        ///number of facts in the program
        uint32_t numFacts = 0;

        ///number of rules in the program
        uint32_t numRules = 0;

        ///rules
        std::list<std::pair<std::list<lit_type>, std::list<lit_type>>> rules;

        ///used for matching literal ids before and after pre processing
        std::unordered_map<uint32_t, uint32_t> atomIdsNewOld;

        ///used for matching literal ids before and after pre processing
        std::unordered_map<uint32_t, uint32_t> atomIdsOldNew;

        ///number of not horn clauses in the program
        uint32_t numNonHornClauses = 0;

        ///number of dual horn clauses in the program
        uint32_t numNonDualHornClauses = 0;

        ///symbol table to match the atom ids with their symbols
        std::unordered_map<uint32_t, const char *> *sTable;

        ///maximal clause size
        uint64_t maxClauseSize = 0;

        ///maximal positive clause size
        uint64_t maxClauseSizePositive = 0;

        ///maximal negative clause size
        uint32_t maxClauseSizeNegative = 0;

        ///occurences of an atom
        std::unordered_map<uint32_t, uint32_t> atomOccurences;

        ///positive occurences of an atom
        std::unordered_map<uint32_t, uint32_t> atomOccurencesPositive;

        ///negative occurences of an atom
        std::unordered_map<uint32_t, uint32_t> atomOccurencesNegative;

        ///variable occurs negative
        std::unordered_map<uint32_t, bool> variableNegative;

        ///variable occurs positive
        std::unordered_map<uint32_t, bool> variablePositive;

        ///negative variable occurrences without helpers
        std::unordered_map<uint32_t, bool> variableNegativeWithoutHelper;

        ///positive variable occurrences without helpers
        std::unordered_map<uint32_t, bool> variablePositiveWithoutHelper;

        ///number of constraints in the program
        uint32_t numConstraints = 0;

        ///maximum size of the positive rules in constraints
        uint32_t maxPositiveRuleSizeConstraint = 0;

        ///maximum size of the positive rules in non constraints
        uint32_t maxPositiveRuleSizeNonConstraint = 0;

        ///num of atom occurences in constraints
        uint32_t atomOccurencesConstraint = 0;

        ///num of atom occurences in non constraints
        uint32_t atomOccurencesNonConstraint = 0;

        ///calculator for the graph Statistics
        GraphStatsCalculator graphStatsCalculator;

        /*
         * private constructors for singleton
         */
        StatsCalculator() : graphStatsCalculator(atomIdsNewOld)
        {
        };

    private:

        /*
         *
         */
        void parseVariableLiteral(std::list<lit_type> body, std::list<lit_type> head);

        /*
         *
         */
        void countAtomOccurences(std::list<lit_type> body, std::list<lit_type> head);

        /*
         * checks if a rule is a constraint
         */
        bool isConstraint(std::list<lit_type> head, std::list<lit_type> body);

        /*
         * calculates the variable occurrences with and without helpers
         */
        void calculateVariables();
    };
}
#endif //CLASP_EXTENDED_STATS_CALCULATOR_H
