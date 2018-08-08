/*
 * FilterableGraph.h
 *
 *  Created on: Jan 30, 2014
 *      Author: Mingxun Wang
 */

#ifndef FILTERABLEGRAPH_H_
#define FILTERABLEGRAPH_H_

#include <set>
#include <map>
#include <vector>
#include <stdlib.h>

#include "Node.h"
#include "Edge.h"
#include "Logger.h"
#include "tuple.h"
#include "BaseGraph.h"

using namespace std;


namespace specnets
{

	class BaseGraph;

        class FilterableGraph: public BaseGraph {
            public:
        
                FilterableGraph();
        
                /**
                * Returns the number of connected components
                */
                int countComponents();
                    
                /**
                * Returns the Edge connecting two nodes, assuming not a multigraph
                */
                Edge * getEdge(Node * from, Node* to);
        
                int getAllComponents(vector<int> & component_sizes, 
                                     vector<unsigned int> & components_node,
                                     vector<vector<Node*> > & nodeLists,
                                     vector<vector<Edge*> > & edgeLists);
                                                      
                int getConnectedComponent(Node * source, vector<Node*> & node_list, vector<Edge*> & edge_list);
                
                virtual int getComponentSizes(vector<int> & component_sizes, vector<unsigned int> & components_node);
            
        };
}

#endif /* FILTERABLEGRAPH_H_ */