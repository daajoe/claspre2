#include <exst/extended_stats_calculator.h>

namespace exst
{
    bool StatsCalculator::isConstraint(std::list<lit_type> head, std::list<lit_type> body)
    {
        std::list<exst::lit_type>::iterator it;
        for (it = head.begin(); it != head.end(); it++)
        {
            if (atomIdsOldNew.count(it->id))
            {
                return true;
            }
        }
        for (it = body.begin(); it != body.end(); it++)
        {
            if (atomIdsOldNew.count(it->id))
            {
                return true;
            }
        }
        return false;
    }

    void StatsCalculator::parseRule(std::list<lit_type> body, std::list<lit_type> head)
    {
        graphStatsCalculator.dependencyGraphStats.addRuleDependencyGraph(body, head);
        graphStatsCalculator.incidenceGraphStats.addRuleIncidenceGraph(body, head);
        std::pair<std::list<lit_type>, std::list<lit_type>> rule;
        rule.first = head;
        rule.second = body;
        rules.push_back(rule);
    }

    void StatsCalculator::calculateStats()
    {
        std::list<std::pair<std::list<exst::lit_type>, std::list<exst::lit_type>>>::iterator it;
        for (it = rules.begin(); it != rules.end(); it++)
        {
            std::list<lit_type> &head = it->first;
            std::list<lit_type> &body = it->second;
            //rules
            numRules++;

            // num facts
            if (body.size() == 0)
                numFacts++;

            uint32_t negative = 0, positive = 0;
            std::list<exst::lit_type>::iterator i;
            for (i = body.begin(); i != body.end(); i++)
            {
                if ((*i).s == -1)
                {
                    negative++;
                } else
                {
                    positive++;
                }
            };
            for (i = head.begin(); i != head.end(); i++)
            {
                if ((*i).s == -1)
                {
                    negative++;
                } else
                {
                    positive++;
                }
            };

            if (isConstraint(head, body))
            {
                maxPositiveRuleSizeConstraint = (uint32_t) (maxPositiveRuleSizeConstraint < positive ? positive
                                                                                                     : maxPositiveRuleSizeConstraint);
                numConstraints++;
                atomOccurencesConstraint += head.size();
                atomOccurencesConstraint += body.size();
            } else
            {
                maxPositiveRuleSizeNonConstraint = (uint32_t) (maxPositiveRuleSizeNonConstraint < positive ? positive
                                                                                                           : maxPositiveRuleSizeNonConstraint);
                atomOccurencesNonConstraint += head.size();
                atomOccurencesNonConstraint += body.size();
            }
            //non horn clause
            if (!((head.size() == 0 && negative <= 1) || (head.size() == 1 && negative <= 0)) && (body.size() > 0))
                ++numNonHornClauses;

            //non dual horn clause
            if ((body.size() - negative > 1) && (body.size() > 0))
                ++numNonDualHornClauses;

            //clause size
            maxClauseSize = (maxClauseSize < body.size() + head.size() ? body.size() + head.size() : maxClauseSize);
            //positive clause size
            maxClauseSizePositive = (maxClauseSizePositive < body.size() - negative + head.size() ? body.size() -
                                                                                                    negative +
                                                                                                    head.size()
                                                                                                  : maxClauseSizePositive);
            //negative clause size
            maxClauseSizeNegative = maxClauseSizeNegative < negative ? negative : maxClauseSizeNegative;

            //atom occurrences
            countAtomOccurences(body, head);

            //number of variables that occur in positive/negative literals
            parseVariableLiteral(body, head);
        }
        calculateVariables();
    }

    void StatsCalculator::countAtomOccurences(std::list<lit_type> body, std::list<lit_type> head)
    {
        std::list<exst::lit_type>::iterator it;
        for (it = body.begin(); it != body.end(); it++)
        {
            uint32_t id = it->id;
            if (it->s == -1)
            {
                atomOccurencesNegative[id]++;
            } else
            {
                atomOccurencesPositive[id]++;
            }
            atomOccurences[id]++;
        }
        for (it = head.begin(); it != head.end(); it++)
        {
            uint32_t id = it->id;
            if (it->s == -1)
            {
                atomOccurencesNegative[id]++;
            } else
            {
                atomOccurencesPositive[id]++;
            }
            atomOccurences[id]++;
        }
    }

    void StatsCalculator::addId(uint32_t before, uint32_t after)
    {
        atomIdsNewOld[after] = before;
        atomIdsOldNew[before] = after;
    }

    void StatsCalculator::printExtendedStats()
    {
        calculateStats();
        std::cout << " \"Extended Stats\" : [\n";

        //dependency graph
        std::cout << "  [\"Dependency Graph Nodes\" , "
                  << graphStatsCalculator.dependencyGraphStats.getDependencyGraph().size() <<
                  "],\n";
        std::cout << "  [\"Dependency Graph Edges\" , "
                  << edgeCount(graphStatsCalculator.dependencyGraphStats.getDependencyGraph()) <<
                  "],\n";

        //incidence graph
        std::cout << "  [\"Incidence Graph Nodes\" , " <<
                  graphStatsCalculator.incidenceGraphStats.getIncidenceGraph().size() << "],\n";
        std::cout << "  [\"Incidence Graph Edges\" , " <<
                  edgeCount(graphStatsCalculator.incidenceGraphStats.getIncidenceGraph()) << "],\n";

        //number of non horn clauses
        std::cout << "  [\"Non Horn Clauses\" , " << numNonHornClauses << "],\n";

        //number of non dual horn clauses
        std::cout << "  [\"Non Dual Horn Clauses\" , " << numNonDualHornClauses << "],\n";

        //maximum clause size
        std::cout << "  [\"max clause size\" , " << maxClauseSize << "],\n";

        //maximum positive/negative clause size, i.e., only positive/negative literals are counted
        std::cout << "  [\"max positive clause size\" , " << maxClauseSizePositive << "],\n";
        std::cout << "  [\"max negative clause size\" , " << maxClauseSizeNegative << "],\n";

        //number of variables that occur as positive/negative literals
        std::cout << "  [\"number of variables that occur as positive literals with helpers\" , " <<
                  variablePositive.size() << "],\n";
        std::cout << "  [\"number of variables that occur as positive literals without helpers\" , " <<
                  variablePositiveWithoutHelper.size() << "],\n";
        std::cout << "  [\"number of variables that occur as negative literals with helpers\" , " <<
                  variableNegative.size() << "],\n";
        std::cout << "  [\"number of variables that occur as negative literals without helpers\" , " <<
                  variableNegativeWithoutHelper.size() << "],\n";

        //maximum positive rule size (constraint/non-constraint)
        std::cout << "  [\"maximum positive rule size constraint\" , " << maxPositiveRuleSizeConstraint << "],\n";
        std::cout << "  [\"maximum positive rule size non-constraint\" , " << maxPositiveRuleSizeNonConstraint <<
                  "],\n";

        //total number of atom occurrences in the program (constraints/non-constraint)
        std::cout << "  [\"total number of atom occurrences constraint\" , " << atomOccurencesConstraint << "],\n";
        std::cout << "  [\"total number of atom occurrences non-constraint\" , " << atomOccurencesNonConstraint <<
                  "],\n";

        //maximum number of occurrences of an atom
        std::cout << "  [\"maximum number of occurrences of an atom\" , " << maxValue(atomOccurences) << "],\n";

        //maximum number of positive occurrences of an atom
        std::cout << "  [\"maximum number of positive occurrences of an atom\" , " <<
                  maxValue(atomOccurencesPositive) << "],\n";

        //maximum number of negative occurrences of an atom
        std::cout << "  [\"maximum number of negative occurrences of an atom\" , " <<
                  maxValue(atomOccurencesNegative) << "]";

        std::list<float>::iterator it;
        for (it = graphStatsCalculator.incidenceGraphStats.getWidths()->begin();
             it != graphStatsCalculator.incidenceGraphStats.getWidths()->end(); it++)
        {
            std::cout << ",\n  [\"reduct width\" , " << (*it) << "]";
        }

        std::cout << "\n ]\n";
        std::flush(std::cout);
    }

    void StatsCalculator::parseVariableLiteral(std::list<lit_type> body, std::list<lit_type> head)
    {
        std::list<exst::lit_type>::iterator it;
        for (it = body.begin(); it != body.end(); it++)
        {
            uint32_t id = it->id;
            if (atomIdsOldNew.count(id) > 0)
            {
                if (it->s == -1)
                {
                    variableNegative[id] = true;
                    variableNegativeWithoutHelper[id] = true;
                } else
                {
                    variablePositive[id] = true;
                    variablePositiveWithoutHelper[id] = true;
                }
            }
        }
        for (it = head.begin(); it != head.end(); it++)
        {
            uint32_t id = it->id;
            if (atomIdsOldNew.count(id) > 0)
            {
                if (it->s == -1)
                {
                    variableNegative[id] = true;
                    variableNegativeWithoutHelper[id] = true;
                } else
                {
                    variablePositive[id] = true;
                    variablePositiveWithoutHelper[id] = true;
                }
            }
        }
    }

    void StatsCalculator::setSymbolTable(std::unordered_map<uint32_t, const char *> &table)
    {
        this->sTable = &table;
    }

    void StatsCalculator::calculateVariables()
    {
        std::list<unsigned int> rem;
        std::unordered_map<unsigned int, bool>::iterator it;
        for (it = variablePositiveWithoutHelper.begin(); it != variablePositiveWithoutHelper.end(); it++)
        {
            if (sTable->count(it->first) == 0)
            {
                rem.push_back((*it).first);
            }
        }
        while (rem.size() > 0)
        {
            variablePositiveWithoutHelper.erase(rem.front());
            rem.pop_front();
        }
        for (it = variableNegativeWithoutHelper.begin(); it != variableNegativeWithoutHelper.end(); it++)
        {
            if (sTable->count((*it).first) == 0)
            {
                rem.push_back((*it).first);
            }
        }
        while (rem.size() > 0)
        {
            variableNegativeWithoutHelper.erase(rem.front());
            rem.pop_front();
        }
    }
}