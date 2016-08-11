#include <exst/incidence_graph_stats.h>

namespace exst
{
    void IncidenceGraphStats::addRuleIncidenceGraph(std::list<lit_type> body, std::list<lit_type> heads)
    {
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, EdgeType>> &myGraph = incidenceGraph;
        std::unordered_map<uint32_t, uint32_t> &avmap = atomVertexMap;
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, EdgeType>> &rbmap = ruleBodyMap;

        uint32_t rule_Vertex = iGraph->addVertex();
        myGraph[rule_Vertex];

        // add body atoms
        std::list<lit_type>::iterator it;
        for (it = body.begin(); it != body.end(); it++)
        {
            uint32_t dId = it->id;
            //if atom is new, add vertex to graph
            if (avmap.count(dId) == 0)
            {
                avmap[dId] = iGraph->addVertex();
                myGraph[avmap[dId]];
            }
            bool neg = it->s == NEGATIVE;
            rbmap[rule_Vertex][avmap[dId]] = neg ? NEG : POS;
            myGraph[avmap[dId]][rule_Vertex] = neg ? NEG : POS;
            myGraph[rule_Vertex][avmap[dId]] = neg ? NEG : POS;
            iGraph->addEdge(avmap[dId], rule_Vertex);
        }

        // add head atoms
        std::list<lit_type>::iterator it_;
        for (it_ = heads.begin(); it_ != heads.end(); it_++)
        {
            uint32_t hId = it_->id;
            //if atom is new, add vertex to graph
            if (avmap.count(hId) == 0)
            {
                avmap[hId] = iGraph->addVertex();
                myGraph[avmap[hId]];
            }
            myGraph[avmap[hId]][rule_Vertex] = HEAD;
            myGraph[rule_Vertex][avmap[hId]] = HEAD;
            iGraph->addEdge(avmap[hId], rule_Vertex);
        }
    }

    void IncidenceGraphStats::resetAssignment()
    {
        if (iGraphReduct == nullptr || iGraph == nullptr)
        {
            return;
        }
        //reset htd reduct graph
        iGraphReduct = iGraph->clone();
        //selectedAtoms.clear();
        ruleBodyMapReduct.insert(ruleBodyMap.begin(), ruleBodyMap.end());
        //reset reduct graph
        incidenceGraphReduct = copyMyGraph(incidenceGraph);
    }

    void IncidenceGraphStats::reduceGraph(lit_type lit)
    {
        if (iGraphReduct == nullptr || iGraph == nullptr)
        {
            return;
        }
        uint32_t nodeIdBody = atomVertexMap[atomIds[lit.id]];
        MyGraph &igraph = incidenceGraphReduct;
        //get all bodies the literal is in
        std::unordered_map<uint32_t, EdgeType> edges = igraph[nodeIdBody];

        std::unordered_map<unsigned int, EdgeType>::iterator it;
        //iterate over all bodies the literal is in
        for (it = edges.begin(); it != edges.end(); it++)
        {
            //if sign is different to the sign of the literal in the body, reduce graph
            if ((*it).second == (lit.s == POSITIVE ? NEG : POS))
            {
                std::unordered_map<uint32_t, EdgeType> bodies = ruleBodyMapReduct[(*it).first];
                std::unordered_map<unsigned int, EdgeType>::iterator bodyIt;
                //remove all edges of wrong body
                for (bodyIt = bodies.begin(); bodyIt != bodies.end(); bodyIt++)
                {
                    //remove only body edges
                    if (bodyIt->second == POS || bodyIt->second == NEG)
                    {
                        igraph[bodyIt->first].erase(it->first);
                        igraph[it->first].erase(bodyIt->first);
                        ruleBodyMapReduct[it->first].erase(bodyIt->first);
                        iGraphReduct->removeEdge(*iGraphReduct->associatedEdgeIds(bodyIt->first, it->first).begin());
                    }
                }
            }
        }
        float w = getTreewidth(iGraphReduct, libraryInstance);
        //calculate and save the Treewidth of the reduct graph
        widths.push_back(w);
    }

    void IncidenceGraphStats::printIGraphReduct()
    {
        if (iGraphReduct == nullptr)
        {
            return;
        }
        std::cout << "{\"_Incidence Graph Reduct_\": [ \n  [\"Nodes\", ";
        std::cout << iGraphReduct->vertexCount();
        std::cout << "],\n  [\"Edges\", ";
        std::cout << iGraphReduct->edgeCount();
        std::cout << "],\n  [\"Treewidth\", ";
        std::cout << widths.back();
        std::cout << "]\n}\n";
    }

    void IncidenceGraphStats::updateAssignment(std::list<lit_type> new_assignment)
    {
        const std::list<lit_type>::iterator &pos = std::search(new_assignment.begin(), new_assignment.end(),
                                                               current_assignment.begin(),
                                                               current_assignment.end());

        if (pos == new_assignment.end())
        {
            //current assignment is not contained in new assignment
            resetAssignment();
            std::list<lit_type>::iterator it;
            for (it = new_assignment.begin(); it != new_assignment.end(); it++)
            {
                reduceGraph(*it);
            }
        } else
        {

            if (new_assignment.size() == current_assignment.size())
            {
                //current assignment = new assignment
                return;
            } else
            {
                //current assignment is contained in new assignment
                std::list<lit_type>::iterator it = new_assignment.begin();
                for (int i = 0; i < current_assignment.size(); i++)
                {
                    ++it;
                }
                for (; it != new_assignment.end(); it++)
                {
                    reduceGraph(*it);
                }
            }
        }
        current_assignment = new_assignment;
    }

    size_t getTreewidth(htd::IMutableMultiGraph *graph, htd::LibraryInstance *libraryInstance)
    {
        if (graph == nullptr || libraryInstance == nullptr)
        {
            return 0;
        }
        //initialize algorithm for decomposition
        htd::IterativeImprovementTreeDecompositionAlgorithm algorithm(libraryInstance,
                                                                      libraryInstance->treeDecompositionAlgorithmFactory().getTreeDecompositionAlgorithm(
                                                                              libraryInstance),
                                                                      WidthMinimizingFitnessFunction());
        // compute decomposition of graph
        htd::ITreeDecomposition *decomposition = algorithm.computeDecomposition(*graph);
        // if there is no decomposition return 0
        return decomposition != nullptr ? decomposition->maximumBagSize() : 0;
    }
}