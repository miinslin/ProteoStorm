/*
 * SpecnetsGraph.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Mingxun Wang
 */


#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "FilterableGraph.h"
#include "Logger.h"
#include "SpecnetsGraph.h"
#include "utils.h"


const int FILTER_REMOVAL_FRACTION = 20;

using namespace std;


namespace specnets
{

  bool EdgeWeightSorter(Edge * left, Edge * right)
  {
    if (left->getWeight() < right->getWeight()) {
      return true;
    }
    return false;
  }

  
  SpecnetsGraph::SpecnetsGraph()
      : FilterableGraph()
  {
  }
  
  SpecnetsGraph::SpecnetsGraph(SpectrumPairSet pairs_set)
      : FilterableGraph()
  {
    for(int i = 0; i < pairs_set.size(); i++) {
      SpectrumPair pair = pairs_set[i];
      int spec1 = pair.spec1;
      int spec2 = pair.spec2;
      
      map<unsigned int, unsigned int>::iterator it;
      
      it=specnets_to_graph_nodes_map.find(spec1);
      
      if(it == specnets_to_graph_nodes_map.end()) {
          Node * new_node = addNode();
          specnets_to_graph_nodes_map[spec1] = new_node->getIndex();
          graph_to_specnets_nodes_map[new_node->getIndex()] = spec1;
      }
      
      it=specnets_to_graph_nodes_map.find(spec2);
      
      if(it == specnets_to_graph_nodes_map.end()) {
          Node * new_node = addNode();
          specnets_to_graph_nodes_map[spec2] = new_node->getIndex();
          graph_to_specnets_nodes_map[new_node->getIndex()] = spec2;
      }
      
      //Add Edge
      Edge * new_edge = addEdge(specnets_to_graph_nodes_map[spec1], specnets_to_graph_nodes_map[spec2]);
      specnets_to_graph_edge_map[i] = new_edge->getIndex();
      graph_to_specnets_edge_map[new_edge->getIndex()] = i;
    }
    graph_pairs_set = pairs_set;
  }
  
  //Filters out edges in connected component by ever increasing score1 to reduce its size to below maximum
  int SpecnetsGraph::filter_graph_component_size(int maximum_component_size, float delta)
  {
     
    bool oversized_component = false;
    
    do {
      oversized_component = false;
        
      vector<int> component_sizes;
      vector<unsigned int> components_node;
      
      //getComponentSizes(component_sizes, components_node);
      vector<vector<Node*> > nodeLists;
      vector<vector<Edge*> > edgeLists;
      getAllComponents(component_sizes, components_node, nodeLists, edgeLists);
      
      for (int i = 0; i < component_sizes.size(); i++) {
        
        if (component_sizes[i] > maximum_component_size) {
          
          DEBUG_MSG("Component [" << i << "] too large. Size = "<< component_sizes[i]);
          oversized_component = true;
          clear_traversal_datastructures();
          
          // Sort the edges by weight and remove a fraction of them (have to remove at least 1)
          std::sort(edgeLists[i].begin(), edgeLists[i].end(), EdgeWeightSorter);
          int removeNumber = max(1, component_sizes[i] / FILTER_REMOVAL_FRACTION);
          DEBUG_VAR(removeNumber);
          
          for (int j = 0; j < removeNumber; j++) {
            //DEBUG_VAR(edgeLists[i][j]->getIndex());
            //Mark this index to be deleted from pairs
            specnets_pairs_deleted.push_back(graph_to_specnets_edge_map[edgeLists[i][j]->getIndex()]);
            removeEdge(edgeLists[i][j]->getIndex());
          }
        }
      }
        
    } while(oversized_component == true);
    
    return 0;
  }
  
  int SpecnetsGraph::get_node_to_connected_component_mapping(map<unsigned int, unsigned int> & mapping){
      
      clear_traversal_datastructures();
      
      unsigned long num_nodes = numNodes();
      
      unsigned int connected_component_count = 0;
      
      for(unsigned int i = 0; i < num_nodes; i++) {
          Node * node1 = getNode(i);
          
          std::map<unsigned int, unsigned int>::iterator it;
          it=mapping.find(graph_to_specnets_nodes_map[node1->getIndex()]);
          if (it != mapping.end()) {
              continue;
          }
          
          vector<Node*> node_list;
          vector<Edge*> edge_list;
          
          connected_component_count++;
          
          if(getConnectedComponent(node1, node_list, edge_list) != -1) {
              for(int j = 0; j < node_list.size(); j++){
                  mapping[graph_to_specnets_nodes_map[node_list[j]->getIndex()]] = connected_component_count;
              }
          }
          else{
              continue;
          }
          
      }
      
      return 0;
  }
  
  float SpecnetsGraph::lowest_score_edge(vector<Edge*> edge_list) {
      float minimum = 1000.f;
      
      for(int i = 0; i < edge_list.size(); i++) {
          //cout<<"edge ptr "<<edge_list[i]<<endl;
          //cout<<"graph edge index "<<edge_list[i]->getIndex()<<endl;
          //cout<<"specnets index "<<graph_to_specnets_edge_map[edge_list[i]->getIndex()]<<endl;
          float score1 = graph_pairs_set[graph_to_specnets_edge_map[edge_list[i]->getIndex()]].score1;
          
          //cout<<i<<"\t"<<score1<<endl;
          
          minimum = min(minimum, score1);
      }
      
      return minimum;
  }
  
  int SpecnetsGraph::remove_low_scoring_edges(vector<Edge*> edge_list, float minimum_threshold){
      
      int number_removed_edges = 0;
      
      for(int i = 0; i < edge_list.size(); i++) {
          //cout<<"edge ptr "<<edge_list[i]<<endl;
          //cout<<"graph edge index "<<edge_list[i]->getIndex()<<endl;
          //cout<<"specnets index "<<graph_to_specnets_edge_map[edge_list[i]->getIndex()]<<endl;
          float score1 = graph_pairs_set[graph_to_specnets_edge_map[edge_list[i]->getIndex()]].score1;
          
          if (score1 < minimum_threshold) {
              //Mark this index to be deleted from pairs
              specnets_pairs_deleted.push_back(graph_to_specnets_edge_map[edge_list[i]->getIndex()]);
              
              //Delete Edge From Graph
              removeEdge(edge_list[i]->getIndex());
              number_removed_edges++;
          }
      }
      
      return number_removed_edges;
  }
  
  vector<unsigned int> SpecnetsGraph::get_pairs_deleted() {
      return specnets_pairs_deleted;
  }
    
}