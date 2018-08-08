/*
 * FilterableGraph.cpp
 *
 *  Created on: Jan 30, 2014
 *      Author: Mingxun Wang
 */


#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "FilterableGraph.h"
#include "utils.h"

using namespace std;


namespace specnets
{
    FilterableGraph::FilterableGraph()
    : BaseGraph()
    {
    }
    
    int FilterableGraph::countComponents(){
        return 0;
    }
    
    int FilterableGraph::getAllComponents(vector<int> & component_sizes, 
                                          vector<unsigned int> & components_node,
                                          vector<vector<Node*> > & nodeLists,
                                          vector<vector<Edge*> > & edgeLists){
        clear_traversal_datastructures();
        
        unsigned long num_nodes = numNodes();
        
        for(unsigned int i = 0; i < num_nodes; i++){
            //TODO yell at adrian to not exit progrma when recoverable
            Node * node1 = getNode(i);
            
            vector<Node*> node_list;
            vector<Edge*> edge_list;
            
            if(getConnectedComponent(node1, node_list, edge_list) != -1){
                component_sizes.push_back(node_list.size());
                components_node.push_back(node1->getIndex());
                nodeLists.push_back(node_list);
                edgeLists.push_back(edge_list);
            }
            else{
                continue;
            }
            
        }
        
        return 0;
    }
    
    int FilterableGraph::getComponentSizes(vector<int> & component_sizes, vector<unsigned int> & components_node){
        clear_traversal_datastructures();
        
        unsigned long num_nodes = numNodes();
        
        for(unsigned int i = 0; i < num_nodes; i++){
            //TODO yell at adrian to not exit progrma when recoverable
            Node * node1 = getNode(i);
            
            vector<Node*> node_list;
            vector<Edge*> edge_list;
            
            if(getConnectedComponent(node1, node_list, edge_list) != -1){
                component_sizes.push_back(node_list.size());
                components_node.push_back(node1->getIndex());
            }
            else{
                continue;
            }
            
        }
        
        return 0;
    }
    
    int FilterableGraph::getConnectedComponent( Node * source, 
                                                vector<Node*> & node_list, 
                                                vector<Edge*> & edge_list){
        if(source == NULL)
            return -1;
        
        if(continueBFS_Undirected(source) == 0){
            return -1;
        }
        
        Node * current_node = source;
        Node * previous_node = source;
        Edge * traversed_edge = NULL;
        do{
            node_list.push_back(current_node);
            //cout<<"NODE " << current_node->getIndex()<<endl;
            //Skip first node
            if(current_node == previous_node){
                current_node = this->nextBFSDFS_Undirected(traversed_edge);
                continue;
            }
            
            
            edge_list.push_back(traversed_edge);
            traversed_edge = NULL;
            
            previous_node = current_node;
            current_node = this->nextBFSDFS_Undirected(traversed_edge);
            
        }while(current_node != NULL);
        
        return 0;
    }
    
    Edge * FilterableGraph::getEdge(Node * from, Node* to){
        if(from == NULL){
            return NULL;
        }
        
        if(to == NULL){
            return NULL;
        }
        
        IndexVector edges = this->getOutEdges(from);
        
        for(int i = 0; i < edges.size(); i++){
            Edge * possible_edge = ((BaseGraph*) this)->getEdge(edges[i]);
            long to_node_index = possible_edge->getToNodeIndex();
            if(to_node_index == to->getIndex()){
                return possible_edge;
            }
        }
        return NULL;
    }
    
    

}