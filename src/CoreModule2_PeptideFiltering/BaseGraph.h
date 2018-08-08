/*
 * BaseGraph.h
 *
 *  Created on: Feb 8, 2012
 *      Author: aguthals
 */

#ifndef BASEGRAPH_H_
#define BASEGRAPH_H_

#include <set>
#include <map>
#include <vector>
#include <stdlib.h>

#include "Node.h"
#include "Edge.h"
#include "Logger.h"
#include "tuple.h"

using namespace std;

namespace specnets
{
  class Node;
  class Edge;
  class EMExcep;
  class ENFExcep;
  class ENFExcep;

  // A Path is defined as a sequence of Edges
  typedef vector<Edge*> Path;

  // A Tree is defined as a collection of Nodes with at most one Path between each pair
  // Each node references the score of the path to a source along with the Edge that points to the Node
  typedef map<Node*, pair<double, Edge*> > Tree;

  typedef vector<Node*> NodeSet;
  typedef vector<Edge*> EdgeSet;

  // A collection of Node or Edge indices
  typedef vector<long> IndexVector;

  // All Node/Edge indices that are vacant
  typedef list<unsigned long> FreeList;

  // References each node to outgoing or incoming edges
  typedef vector<pair<IndexVector, FreeList> > AdjacList;

  class BaseGraph
  {
  public:

    // binary format version
    static const unsigned short BIN_VERSION;

    // binary format sub-version
    static const unsigned short BIN_SUBVERSION;

    // string identifier for version
    static const string BIN_VERSION_ID;

    // string identifier for sub-version
    static const string BIN_SUBVERSION_ID;

    /**
     * Default constructor. Although no nodes or edges are constructed, this
     *   allocates memory for an initial number of nodes and edges
     * @param numNodes number of nodes to allocate memory for
     * @param numEdges number of edges to allocate memory for
     */
    BaseGraph(unsigned long numNodes = 0, unsigned long numEdges = 0);

    /**
     * Copy constructor
     */
    BaseGraph(const BaseGraph& other);

    /**
     * De-constructor
     */
    virtual ~BaseGraph(void);

    /**
     * Returns a string identifier of this graph type w/ its memory address
     */
    virtual string toString(void) const;

    /**
     * Creates a new Node on the heap as a copy of another. This is called internally
     *   when ever a node is copied.
     * @param copyNode
     * @return pointer to new node
     */
    virtual Node* cloneNode(Node& copyNode) const;

    /**
     * Creates a new Edge on the heap as a copy of another. This is called internally
     *   when ever an edge is copied.
     * @param copyEdge
     * @return pointer to new edge
     */
    virtual Edge* cloneEdge(Edge& copyEdge) const;

    /**
     * Creates a new Node on the heap. This is called internally
     *   when ever a Node is created.
     * @return pointer to new node
     */
    virtual Node* createNode(void) const;

    /**
     * Creates a new Edge on the heap. This is called internally
     *   when ever an Edge is created.
     * @return pointer to new edge
     */
    virtual Edge* createEdge(void) const;

    /**
     * Checks if this graph's internal data structures are
     *   consistent with each other and each Node and Edge data
     *   structures. If they are not consistent, an abort signal is
     *   declared with an error message
     */
    virtual void validateGraph(void) const;

    /**
     * Checks if an Edge's to/from pointers are consistent with each
     *   Node. If they are not consistent, an abort signal is declared
     *   with an error message
     * @param edge
     */
    virtual void validateEdge(Edge* edge) const;

    /**
     * Copies the structure and data of another graph into this graph.
     *   Nodes and Edges are copied using cloneNode and cloneEdge.
     * @param other
     * @return reference to this graph
     */
    BaseGraph &operator=(const BaseGraph &other);

    /**
     * Same as calling operator=(), except no reference is returned.
     */
    virtual void copy(BaseGraph& other);

    /**
     * Appends the entire structure of another graph to this one. Node
     *   and Edge indices from otherGraph are updated to make room.
     * @param otherGraph
     */
    virtual void appendGraph(const BaseGraph& otherGraph);

    /**
     * Saves this graph structure as a text file that is readable by
     *   graphviz.
     * @param filename
     * @return true if file was saved successfully, false if not.
     */
    virtual bool saveGraphviz(const string& filename) const;

    /**
     * Saves the graph structure and all associated data in binary format
     * @param filename
     * @return true if file was saved successfully, false if not.
     */
    virtual bool saveBinaryFile(const string& filename);

    /**
     * Adds the current binary format version and subversion to a map
     */
    virtual void
    addBinaryVersionInfo(map<string, unsigned short>& versions) const;

    /**
     * Loads the graph structure and all associated data in binary format
     * @param filename
     * @return true if file was loaded successfully, false if not.
     */
    virtual bool loadBinaryFile(const string& filename);

    /**
     * Saves the graph structure and all associated data to a binary stream
     *   (called by saveBinaryFile())
     * @param fp pointer to binary file stream
     * @return true if all data was successfully saved, false if not
     */
    virtual bool saveGraphToBinaryStream(FILE* fp);

    /**
     * Loads the graph structure and all associated data from a binary stream
     *   (called by loadBinaryFile())
     * @param fp pointer to binary file stream
     * @param versions identifies the loaded version/subversion of each class
     *   that was in the file (ie all Node and Edge descendant classes)
     * @return true if all data was successfully loaded, false if not
     */
    virtual bool
    loadGraphFromBinaryStream(FILE* fp, map<string, unsigned short>& versions);

    /**
     * Called internally when Nodes/Edges are modified to indicate that the graph
     *   was changed so BFS/DFS states are marked invalid
     */
    virtual void modifiedGraph();

    /**
     * constructor method
     */
    void initialize(unsigned long numNodes = 0, unsigned long numEdges = 0);

    /**
     * Creates a new node by calling cloneNode() with a template node
     * @param copyNode reference to node to copy
     * @return pointer to new node on the heap
     */
    Node* addNode(Node& copyNode);

    /**
     * Creates a new node by calling createNode()
     * @return pointer to new node on the heap
     */
    Node* addNode(void);

    /**
     * Creates a new edge by calling createEdge(), then has
     *   it point from an origin node to a destination node
     * @param originIdx index of orgin Node
     * @param destIdx index of destination Node
     * @return pointer to new node on the heap
     */
    Edge* addEdge(const unsigned long& originIdx, const unsigned long& destIdx);

    /**
     * Creates a new edge by calling cloneEdge() on a template, then has
     *   it point from an origin node to a destination node
     * @param originIdx index of orgin Node
     * @param destIdx index of destination Node
     * @param copyEdge edge to copy to new edge
     * @return pointer to new edge on the heap
     */
    Edge* addEdge(const unsigned long& originIdx,
                  const unsigned long& destIdx,
                  Edge& copyEdge);

    /**
     * Moves an existing edge to a new location
     * @param edgeIdx index of edge to move
     * @param originIdx index of new origin node
     * @param destIdx index of new destination node
     * @return pointer to Edge that was moved
     */
    Edge* moveEdge(const unsigned long& edgeIdx,
                   const unsigned long& originIdx,
                   const unsigned long& destIdx);

    /**
     * Creates a new edge by calling createEdge(), then has
     *   it point from an origin node to a destination node
     * @param origin pointer to orgin Node
     * @param dest pointer to destination Node
     * @return pointer to new edge on the heap
     */
    inline Edge* addEdge(Node* origin, Node* dest)
    {
      return addEdge(origin->getIndex(), dest->getIndex());
    }

    /**
     * Creates a new edge by calling cloneEdge() on a template, then has
     *   it point from an origin node to a destination node
     * @param origin pointer to orgin Node
     * @param dest pointer to destination Node
     * @param copyEdge edge to copy to new edge
     * @return pointer to new edge on the heap
     */
    inline Edge* addEdge(Node* origin, Node* dest, Edge& copyEdge)
    {
      return addEdge(origin->getIndex(), dest->getIndex(), copyEdge);
    }

    /**
     * Removes a node at a given index from this graph (also frees its memory)
     * @param nodeIdx index of Node to delete
     */
    void removeNode(unsigned long nodeIdx);

    /**
     * Removes a node from this graph (also frees its memory)
     * @param node ponter to a pointer of a Node to delete. This will point to
     *   zero upon completion.
     */
    inline void removeNode(Node** node)
    {
      removeNode((*node)->getIndex());
      *node = 0;
    }

    /**
     * Removes an edges at a given index from this graph (also frees its memory)
     * @param edgeIdx index of Edge to delete
     */
    void removeEdge(unsigned long edgeIdx);

    /**
     * Removes an edge from this graph (also frees its memory)
     * @param edge ponter to a pointer of an Edge to delete. This will point to
     *   zero upon completion.
     */
    inline void removeEdge(Edge** edge)
    {
      removeEdge((*edge)->getIndex());
      *edge = 0;
    }

    /**
     * Returns a pointer to the Node at a given index. This will abort if
     *   the index is out-of-bounds
     * @param nodeIndex
     * @return Node pointer
     */
    inline Node* getNode(unsigned long nodeIndex) const
    {
      if (nodeIndex >= m_nodes.size() || m_nodes[nodeIndex] == 0)
      {
        ERROR_MSG("Cannot find node " << nodeIndex);
        abort();
      }
      return m_nodes[nodeIndex];
    }

    /**
     * Returns a pointer to the Node at a given index.
     * @param nodeIndex
     * @return Node pointer
     */
    inline Edge* getEdge(unsigned long edgeIndex) const
    {
      if (edgeIndex >= m_edges.size() || m_edges[edgeIndex] == 0)
      {
        ERROR_MSG("Cannot find edge " << edgeIndex);
        abort();
      }
      return m_edges[edgeIndex];
    }

    /**
     * Same as getNode but with no index-out-of-bounds check
     */
    inline Node* lookupNode(unsigned long nodeIndex) const
    {
      return m_nodes[nodeIndex];
    }

    /**
     * Same as getEdge but with no index-out-of-bounds check
     */
    inline Edge* lookupEdge(unsigned long edgeIndex) const
    {
      return m_edges[edgeIndex];
    }

    /**
     * Checks if a Node is in this graph
     * @param node
     * @return true if node is found in this graph, false if not
     */
    inline bool containsNode(Node* node) const
    {
      return (node != 0) && (node->getIndex() < m_nodes.size())
          && (m_nodes[node->getIndex()] == node);
    }

    /**
     * Checks if an Edge is in this graph
     * @param edge
     * @return true if edge is found in this graph, false if not
     */
    inline bool containsEdge(Edge* edge) const
    {
      return (edge != 0) && (edge->getIndex() < m_edges.size())
          && (m_edges[edge->getIndex()] == edge);
    }

    /**
     * Returns a vector of all edges originating from a node
     * @param nodeIndex index of Node to find edges from
     * @return constant reference to vector of outgoing edges. This is
     *   a vector of edge indices, so some values may be negative and invalid
     */
    inline const IndexVector& getOutEdges(unsigned long nodeIndex) const
    {
      if (nodeIndex >= m_nodes.size() || m_nodes[nodeIndex] == 0)
      {
        ERROR_MSG("Cannot find node " << nodeIndex);
        abort();
      }
      return m_outEdges[nodeIndex].first;
    }

    /**
     * Returns a vector of all edges originating from a node
     * @param node pointer to Node to find edges from
     * @return constant reference to vector of outgoing edges. This is
     *   a vector of edge indices, so some values may be negative and invalid
     */
    inline const IndexVector& getOutEdges(const Node* node) const
    {
      return getOutEdges(node->m_index);
    }

    /**
     * Returns a vector of all edges pointing to a node
     * @param nodeIndex index of Node to find edges to
     * @return constant reference to vector of incoming edges. This is
     *   a vector of edge indices, so some values may be negative and invalid
     */
    inline const IndexVector& getInEdges(unsigned long nodeIndex) const
    {
      if (nodeIndex >= m_nodes.size() || m_nodes[nodeIndex] == 0)
      {
        ERROR_MSG("Cannot find node " << nodeIndex);
        abort();
      }
      return m_inEdges[nodeIndex].first;
    }

    /**
     * Returns a vector of all edges pointing to a node
     * @param node pointer to a Node to find edges to
     * @return constant reference to vector of incoming edges. This is
     *   a vector of edge indices, so some values may be negative and invalid
     */
    inline const IndexVector& getInEdges(const Node* node) const
    {
      return getInEdges(node->m_index);
    }

    /**
     * Returns the number of valid outgoing edges from a given node
     * @param node
     * @return number of valid outgoing edges
     */
    inline unsigned long getNumOutEdges(const Node* node) const
    {
      return getOutEdges(node).size() - m_outEdges[node->m_index].second.size();
    }

    /**
     * Returns the number of valid incoming edges from a given node
     * @param node
     * @return number of valid incoming edges
     */
    inline unsigned long getNumInEdges(const Node* node) const
    {
      return getInEdges(node).size() - m_inEdges[node->m_index].second.size();
    }

    /**
     * Returns the number of valid edges in this graph
     */
    inline unsigned long numEdges(void) const
    {
      return m_edges.size() - m_freeEdges.size();
    }

    /**
     * Returns the number of valid nodes in this graph
     */
    inline unsigned long numNodes(void) const
    {
      return m_nodes.size() - m_freeNodes.size();
    }

    /**
     * Returns the maximum Edge index that is currently available
     */
    inline unsigned long maxEdgeIdx(void) const
    {
      return m_edges.size() - 1;
    }

    /**
     * Returns the maximum Node index that is currently available
     */
    inline unsigned long maxNodeIdx(void) const
    {
      return m_nodes.size() - 1;
    }

    /**
     * Returns the string label of this graph
     */
    inline string getLabel(void) const
    {
      return m_label;
    }

    /**
     * Compresses all node and edge indices
     * @param outputNewNodeIdxs if specified, this will contain a mapping
     *   of old node indices to new indices
     * @param outputNewEdgeIdxs if specified, this will contain a mapping
     *   of old edge indices to new indices
     */
    virtual void compress(vector<long>* outputNewNodeIdxs = 0,
                          vector<long>* outputNewEdgeIdxs = 0);

    /**
     * Returns a tree (originating from source) encoding the
     *   lightest path from source to every other node
     * @param source
     * @param outputPaths
     */
    void getLightestPaths(Node* source, Tree& outputPaths);

    /**
     * Computes the heaviest path from source to sink. The graph must be a
     *   DAG for this work in polynomial time.
     * @param source
     * @param sink
     * @param outputPath
     * @return pair.first-> true of the graph is a DAG and a heaviest path was found
     *         pair.second -> score of the path from source to sink
     */
    pair<bool, double> getHeaviestPathDAG(Node* source,
                                          Node* sink,
                                          Path& outputPath);

    /**
     * Gets a topological ordering of all nodes
     * @return true if this graph is a DAG and an ordering was found.
     */
    bool getTopologicalOrderingDAG(list<Node*>& outputOrder);

    /**
     * Begins a breadth-first search from source
     * @param source
     * @return pointer to the first node in the search
     */
    Node* beginBFS(Node* source);

    /**
     * Begins a depth-first search from source
     * @param source
     * @return pointer to the next first in the search
     */
    Node* beginDFS(Node* source);
    
    /**
     * Begins a breadth-first search from source in Undirected Graph
     * @param source
     * @return pointer to the first node in the search
     */
    Node* beginBFS_Undirected(Node* source);

    /**
     * Begins a depth-first search from source, in Undirected Graph
     * @param source
     * @return pointer to the next first in the search
     */
    Node* beginDFS_Undirected(Node* source);
    
    /**
     * Continues a breadth-first search from source in Undirected Graph
     * @param source
     * @return pointer to the first node in the search
     */
    Node* continueBFS_Undirected(Node* new_source);

    /**
     * Begins a depth-first search from source, in Undirected Graph
     * @param source
     * @return pointer to the next first in the search
     */
    Node* continueDFS_Undirected(Node* new_source);
    
    

    /**
     * Called to get the next discovered node in either
     *   BFS or DFS, If NULL, then there are no more nodes.
     */
    Node* nextBFSDFS(void);
    
    /**
     * Called to get the next discovered node in either
     *   BFS or DFS, treating graph as undirected, 
     *   If NULL, then there are no more nodes.
     * @param traversed_edge
     */
    Node* nextBFSDFS_Undirected(Edge * &traversed_edge);

  protected:
    string m_label;
    
    void clear_traversal_datastructures();

  private:

    NodeSet m_nodes;
    EdgeSet m_edges;

    FreeList m_freeNodes;
    FreeList m_freeEdges;

    AdjacList m_outEdges;
    AdjacList m_inEdges;

    list<Edge*> m_queueStack;
    list<int> m_queueStack_undirected_traversedirection; //0 is from, 1 is to, describing originating node
    set<Node*> m_visited;
    bool m_QueueIsStack;
    bool m_okTraverse;

    void m_beginTraverse(Node* source, bool undirected = false);
    void m_continueTraverse(Node* source, bool undirected = false);

    void m_clearNodes(void);

    unsigned long m_getFreeNodeIdx(void);

    unsigned long m_getFreeEdgeIdx(void);

    unsigned long m_getFreeInAdacIdx(unsigned long nodeIdx);

    unsigned long m_getFreeOutAdacIdx(unsigned long nodeIdx);

  };
}

#endif /* BASEGRAPH_H_ */
