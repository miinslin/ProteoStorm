/*
 * SpecnetsGraph.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Mingxun Wang
 */

#ifndef SPECNETSGRAPH_H_
#define SPECNETSGRAPH_H_

#include <set>
#include <map>
#include <vector>
#include <stdlib.h>

#include "Node.h"
#include "Edge.h"
#include "Logger.h"
#include "tuple.h"
#include "FilterableGraph.h"
#include "SpectrumPairSet.h"

using namespace std;


namespace specnets
{

        class BaseGraph;
        class FilterableGraph;

        class SpecnetsGraph: public FilterableGraph {
            public:
                SpecnetsGraph();
                
                SpecnetsGraph(SpectrumPairSet pairs_set);
                
                int filter_graph_component_size(int maximum_component_size, float delta = 0.1);
                
                int get_node_to_connected_component_mapping(map<unsigned int, unsigned int> & mapping);
                
                vector<unsigned int> get_pairs_deleted();
                
            private:
                map<unsigned int, unsigned int> specnets_to_graph_nodes_map;
                map<unsigned int, unsigned int> specnets_to_graph_edge_map;
                
                map<unsigned int, unsigned int> graph_to_specnets_nodes_map;
                map<unsigned int, unsigned int> graph_to_specnets_edge_map;
                
                SpectrumPairSet graph_pairs_set;
                
                vector<unsigned int> specnets_pairs_deleted;
                
                float lowest_score_edge(vector<Edge*> edge_list);
                
                /**
                 * Removes Edges below minimum threshold and return number of edges removed. 
                **/
                int remove_low_scoring_edges(vector<Edge*> edge_list, float minimum_threshold);
            
        };
}

#endif /* SPECNETSGRAPH_H_ */