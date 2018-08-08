/*
 * BaseGraph.cpp
 *
 *  Created on: Feb 9, 2012
 *      Author: aguthals
 */

#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "HashHeap.h"
#include "BaseGraph.h"
#include "utils.h"

using namespace std;

namespace specnets
{

  const unsigned short BaseGraph::BIN_VERSION = 1;
  const unsigned short BaseGraph::BIN_SUBVERSION = 1;

  const string BaseGraph::BIN_VERSION_ID = "BaseGraph_binVersion";
  const string BaseGraph::BIN_SUBVERSION_ID = "BaseGraph_binSubVersion";

  BaseGraph::BaseGraph(unsigned long numNodes, unsigned long numEdges) :
    m_label("Graph"), m_nodes(numNodes), m_edges(numEdges), m_freeNodes(),
        m_freeEdges(), m_outEdges(numNodes), m_inEdges(numNodes),
        m_queueStack(), m_visited(), m_QueueIsStack(false), m_okTraverse(false)
  {
    this->initialize(numNodes, numEdges);
  }

  BaseGraph::BaseGraph(const BaseGraph& other) :
    m_label("Graph"), m_nodes(other.m_nodes.size()),
        m_edges(other.m_edges.size()), m_freeNodes(), m_freeEdges(),
        m_outEdges(other.m_outEdges.size()), m_inEdges(other.m_inEdges.size()),
        m_queueStack(), m_visited(), m_QueueIsStack(false), m_okTraverse(false)
  {
    this->operator =(other);
  }

  BaseGraph::~BaseGraph(void)
  {
    try
    {
      this->validateGraph();
    }
    catch (exception& e)
    {
      ERROR_MSG("Graph is inconsistent at deallocation! - trying to free memory anyways");
    }
    m_clearNodes();
  }

  Node* BaseGraph::cloneNode(Node& copyNode) const
  {
    return new Node(copyNode);
  }

  Edge* BaseGraph::cloneEdge(Edge& copyEdge) const
  {
    return new Edge(copyEdge);
  }

  Node* BaseGraph::createNode(void) const
  {
    return new Node();
  }

  Edge* BaseGraph::createEdge(void) const
  {
    return new Edge();
  }

  void BaseGraph::validateGraph(void) const
  {

    for (long i = 0; i < m_nodes.size(); i++)
    {
      if (m_nodes[i] == 0)
        continue;

      if (m_nodes[i]->getIndex() != i)
      {
        ERROR_MSG("Invalid node index " << m_nodes[i]->getIndex() << " at index " << i);
        abort();
      }
    }

    for (long i = 0; i < m_edges.size(); i++)
    {
      if (m_edges[i] == 0)
        continue;

      if (m_edges[i]->getIndex() != i)
      {
        ERROR_MSG("Invalid edge index " << m_nodes[i]->getIndex() << " at index " << i);
        abort();
      }
      validateEdge(m_edges[i]);
    }

    for (long i = 0; i < m_nodes.size(); i++)
    {
      if (m_nodes[i] == 0)
        continue;
      for (long j = 0; j < m_outEdges[i].first.size(); j++)
      {
        if (m_outEdges[i].first[j] < 0)
          continue;

        Edge* edge = this->getEdge(m_outEdges[i].first[j]);
        validateEdge(edge);
      }
      for (long j = 0; j < m_inEdges[i].first.size(); j++)
      {
        if (m_inEdges[i].first[j] < 0)
          continue;

        Edge* edge = this->getEdge(m_inEdges[i].first[j]);
        validateEdge(edge);
      }
    }
  }

  void BaseGraph::validateEdge(Edge* edge) const
  {
    if (m_nodes[edge->m_fromNodeIndex] == 0)
    {
      ERROR_MSG("Invalid from node index " << edge->m_fromNodeIndex << " for edge " << edge->m_index);
      abort();
    }
    if (m_nodes[edge->m_toNodeIndex] == 0)
    {
      ERROR_MSG("Invalid to node index " << edge->m_toNodeIndex << " for edge " << edge->m_index);
      abort();
    }

    if (m_outEdges[edge->m_fromNodeIndex].first[edge->m_fromIndex]
        != edge->getIndex())
    {
      ERROR_MSG("Invalid index " << edge->m_fromIndex << " in outgoing edge vector " << edge->m_fromNodeIndex << " for edge " << edge->m_index);
      DEBUG_VAR(m_outEdges[edge->m_fromNodeIndex].first[edge->m_fromIndex]);
      abort();
    }

    if (m_inEdges[edge->m_toNodeIndex].first[edge->m_toIndex]
        != edge->getIndex())
    {
      ERROR_MSG("Invalid index " << edge->m_toIndex << " in incoming edge vector " << edge->m_toNodeIndex << " for edge " << edge->m_index);
      abort();
    }
  }

  string BaseGraph::toString(void) const
  {
    ostringstream out;
    out << m_label << " at <" << this << ">";
    return out.str();
  }

  BaseGraph & BaseGraph::operator=(const BaseGraph &other)
  {
    if (this == &other)
    {
      return *this;
    }

    this->initialize(other.m_nodes.size(), other.m_edges.size());

    m_freeNodes = other.m_freeNodes;
    for (unsigned long i = 0; i < other.m_nodes.size(); i++)
    {
      Node* otherNode = other.m_nodes[i];
      if (otherNode != 0)
      {
        m_nodes[i] = cloneNode(*otherNode);
        m_outEdges[i] = other.m_outEdges[i];
        m_inEdges[i] = other.m_inEdges[i];
      }
    }

    m_freeEdges = other.m_freeEdges;
    for (unsigned long i = 0; i < other.m_edges.size(); i++)
    {
      Edge* otherEdge = other.m_edges[i];
      if (otherEdge != 0)
      {
        m_edges[i] = cloneEdge(*otherEdge);
      }
    }

    return *this;
  }

  void BaseGraph::copy(BaseGraph& other)
  {
    this->operator =(other);
  }

  void BaseGraph::appendGraph(const BaseGraph& other)
  {
    map<unsigned long, unsigned long> indexMap;
    for (unsigned long i = 0; i < other.m_nodes.size(); i++)
    {
      Node* otherNode = other.m_nodes[i];
      if (otherNode != 0)
      {
        Node* myCopy = this->addNode(*otherNode);
        indexMap[i] = myCopy->getIndex();
      }
    }

    for (unsigned long i = 0; i < other.m_edges.size(); i++)
    {
      Edge* otherEdge = other.m_edges[i];
      if (otherEdge != 0)
      {
        unsigned long fromIdx = indexMap[otherEdge->getFromNodeIndex()];
        unsigned long toIdx = indexMap[otherEdge->getToNodeIndex()];
        Edge* newEdge = this->addEdge(fromIdx, toIdx, *otherEdge);
      }
    }
  }

  bool BaseGraph::saveGraphviz(const string& filename) const
  {
    validateGraph();

    ofstream out(filename.c_str(), ios::binary);
    if (!out)
    {
      ERROR_MSG("Error opening \'" << filename << "\' for writing !");
      return false;
    }
    out << "digraph G {\n";

    for (unsigned long i = 0; i < m_nodes.size(); i++)
    {
      Node* nodePtr = m_nodes[i];
      if (nodePtr == 0)
      {
        continue;
      }

      out << nodePtr->getIndex() << "\t[label=\""
          << nodePtr->getGraphvizLabel() << "\", style=filled";

      if (nodePtr->getGraphvizFillColor().length() > 0)
      {
        out << ", fillcolor=" << nodePtr->getGraphvizFillColor();
      }
      out << "];\n";
    }

    for (unsigned long i = 0; i < m_edges.size(); i++)
    {
      Edge* edgePtr = m_edges[i];
      if (edgePtr == 0)
      {
        continue;
      }
      out << edgePtr->getFromNodeIndex() << " -> " << edgePtr->getToNodeIndex()
          << " [label = \"" << edgePtr->getGraphvizLabel() << "\"];\n";
    }
    out << "}";

    out.close();
    return true;
  }

  bool BaseGraph::saveBinaryFile(const string& filename)
  {
    FILE* fp = fopen(filename.c_str(), "wb");
    if (fp == 0)
    {
      ERROR_MSG("Error opening \'" << filename << "\' for writing !");
      return false;
    }

    map<string, unsigned short> versions;
    addBinaryVersionInfo(versions);

    if (!writeStringMapToBinaryStream<unsigned short> (fp, versions))
    {
      ERROR_MSG("Error saving version info for \'" << toString() << "\'");
      fclose(fp);
      return false;
    }

    if (!saveGraphToBinaryStream(fp))
    {
      ERROR_MSG("Error saving \'" << toString() << "\'");
      fclose(fp);
      return false;
    }
    fclose(fp);
    return true;
  }

  void BaseGraph::addBinaryVersionInfo(map<string, unsigned short>& versions) const
  {
    versions[BIN_VERSION_ID] = BIN_VERSION;
    versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
    Node* temp = createNode();
    temp->addBinaryVersionInfo(versions);
    Edge* tempE = createEdge();
    tempE->addBinaryVersionInfo(versions);

    delete temp;
    delete tempE;
  }

  bool BaseGraph::loadBinaryFile(const string& filename)
  {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (fp == 0)
    {
      ERROR_MSG("Error opening \'" << filename << "\' for reading !");
      return false;
    }

    map<string, unsigned short> versions;
    if (!readStringMapFromBinaryStream<unsigned short> (fp, versions))
    {
      ERROR_MSG("Error reading version info for \'" << toString() << "\'");
      fclose(fp);
      return false;
    }

    if (!loadGraphFromBinaryStream(fp, versions))
    {
      ERROR_MSG("Error loading \'" << toString() << "\'");
      fclose(fp);
      return false;
    }
    fclose(fp);
    return true;
  }

  bool BaseGraph::saveGraphToBinaryStream(FILE* fp)
  {
    validateGraph();
    compress();

    unsigned int count;
    unsigned long nNodes = numNodes();
    unsigned long nEdges = numEdges();

    count = fwrite(&nNodes, sizeof(unsigned long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the number of nodes");
      return false;
    }
    count = fwrite(&nEdges, sizeof(unsigned long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the number of edges");
      return false;
    }

    unsigned int labelLength = m_label.length() + 1;
    count = fwrite(&labelLength, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the length of the graph label");
      return false;
    }

    char* labelArr = (char*)malloc(sizeof(char) * labelLength);
    for (unsigned int i = 0; i < labelLength - 1; i++)
    {
      labelArr[i] = m_label.at(i);
    }
    labelArr[labelLength - 1] = '\0';

    count = fwrite(labelArr, sizeof(char), labelLength, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the graph label");
      free(labelArr);
      return false;
    }
    free(labelArr);

    for (unsigned long i = 0; i < m_nodes.size(); i++)
    {
      Node* nodePtr = m_nodes[i];
      if (!nodePtr->saveToBinaryStream(fp))
      {
        ERROR_MSG("Could not write node \'" << nodePtr->toString() << "\'");
        return false;
      }
    }
    for (unsigned long i = 0; i < m_edges.size(); i++)
    {
      Edge* edgePtr = m_edges[i];
      if (!edgePtr->saveToBinaryStream(fp))
      {
        ERROR_MSG("Could not write edge \'" << edgePtr->toString() << "\'");
        return false;
      }
    }

    return true;
  }

  bool BaseGraph::loadGraphFromBinaryStream(FILE* fp,
                                            map<string, unsigned short>& versions)
  {
    if (fp == 0)
    {
      ERROR_MSG("File pointer is zero");
      return false;
    }

    unsigned int count;
    unsigned long nNodes = 0;
    unsigned long nEdges = 0;
    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    count = fread(&nNodes, sizeof(unsigned long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the number of nodes");
      return false;
    }
    count = fread(&nEdges, sizeof(unsigned long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the number of edges");
      return false;
    }

    initialize(nNodes, nEdges);

    unsigned int labelLength = 0;
    count = fread(&labelLength, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could read the length of the graph label");
      return false;
    }

    char* labelArr = (char*)malloc(sizeof(char) * labelLength);

    count = fread(labelArr, sizeof(char), labelLength, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the graph label");
      free(labelArr);
      return false;
    }
    m_label = labelArr;
    free(labelArr);

    for (unsigned int i = 0; i < nNodes; i++)
    {
      Node* newNode = addNode();
      if (!newNode->loadFromBinaryStream(fp, versions))
      {
        ERROR_MSG("Could not read node #" << i);
        return false;
      }
    }

    for (unsigned int i = 0; i < nEdges; i++)
    {
      Edge* newEdge = createEdge();
      if (!newEdge->loadFromBinaryStream(fp, versions))
      {
        ERROR_MSG("Could not read edge #" << i);
        return false;
      }
      Edge* addedEdge = addEdge(newEdge->getFromNodeIndex(),
                                newEdge->getToNodeIndex(),
                                *newEdge);
    }

    return true;
  }

  void BaseGraph::modifiedGraph()
  {
    m_okTraverse = false;
  }

  void BaseGraph::initialize(unsigned long numNodes, unsigned long numEdges)
  {
    modifiedGraph();

    m_clearNodes();

    m_nodes.assign(numNodes, (Node*)0);
    m_edges.assign(numEdges, (Edge*)0);

    pair<IndexVector, FreeList> emptyPair;
    m_outEdges.assign(numNodes, emptyPair);
    m_inEdges.assign(numNodes, emptyPair);

    m_freeNodes.clear();
    for (unsigned long i = 0; i < numNodes; i++)
    {
      m_freeNodes.push_back(i);
    }

    m_freeEdges.clear();
    for (unsigned long i = 0; i < numEdges; i++)
    {
      m_freeEdges.push_back(i);
    }
  }

  Node* BaseGraph::addNode(Node& copyNode)
  {
    unsigned long nodeIdx = this->m_getFreeNodeIdx();
    Node* newNode = cloneNode(copyNode);
    newNode->m_index = nodeIdx;
    m_nodes[nodeIdx] = newNode;
    return newNode;
  }

  Node* BaseGraph::addNode(void)
  {
    unsigned long nodeIdx = this->m_getFreeNodeIdx();
    Node* newNode = createNode();
    newNode->m_index = nodeIdx;
    m_nodes[nodeIdx] = newNode;
    return newNode;
  }

  Edge* BaseGraph::addEdge(const unsigned long& originIdx,
                           const unsigned long& destIdx)
  {
    if (originIdx >= m_nodes.size() || destIdx >= m_nodes.size())
    {
      ERROR_MSG("Invalid node indices " << originIdx << " and " << destIdx);
      abort();
    }
    if (m_nodes[originIdx] == 0 || m_nodes[destIdx] == 0)
    {
      ERROR_MSG("One of " << originIdx << " and " << destIdx << " does not exist as a node");
      abort();
    }
    if (originIdx == destIdx)
    {
      DEBUG_TRACE;
      abort();
    }

    unsigned long edgeIdx = this->m_getFreeEdgeIdx();
    unsigned long inIdx = this->m_getFreeInAdacIdx(destIdx);
    unsigned long outIdx = this->m_getFreeOutAdacIdx(originIdx);

    Edge* newEdge = createEdge();
    newEdge->m_index = edgeIdx;
    newEdge->m_fromIndex = outIdx;
    newEdge->m_toIndex = inIdx;
    newEdge->m_fromNodeIndex = originIdx;
    newEdge->m_toNodeIndex = destIdx;

    m_outEdges[originIdx].first[outIdx] = edgeIdx;
    m_inEdges[destIdx].first[inIdx] = edgeIdx;

    m_edges[edgeIdx] = newEdge;

    return newEdge;
  }

  Edge* BaseGraph::addEdge(const unsigned long& originIdx,
                           const unsigned long& destIdx,
                           Edge& copyEdge)
  {
    if (originIdx >= m_nodes.size() || destIdx >= m_nodes.size())
    {
      ERROR_MSG("Invalid node indices " << originIdx << " and " << destIdx);
      abort();
    }
    if (m_nodes[originIdx] == 0 || m_nodes[destIdx] == 0)
    {
      ERROR_MSG("One of " << originIdx << " and " << destIdx << " does not exist as a node");
      abort();
    }
    if (originIdx == destIdx)
    {
      ERROR_MSG("Origin and destination nodes are the same!!");
      abort();
    }

    unsigned long edgeIdx = this->m_getFreeEdgeIdx();
    unsigned long inIdx = this->m_getFreeInAdacIdx(destIdx);
    unsigned long outIdx = this->m_getFreeOutAdacIdx(originIdx);

    Edge* newEdge = cloneEdge(copyEdge);
    newEdge->m_index = edgeIdx;
    newEdge->m_fromIndex = outIdx;
    newEdge->m_toIndex = inIdx;
    newEdge->m_fromNodeIndex = originIdx;
    newEdge->m_toNodeIndex = destIdx;

    m_outEdges[originIdx].first[outIdx] = edgeIdx;
    m_inEdges[destIdx].first[inIdx] = edgeIdx;

    m_edges[edgeIdx] = newEdge;

    return newEdge;
  }

  Edge* BaseGraph::moveEdge(const unsigned long& edgeIdx,
                            const unsigned long& originIdx,
                            const unsigned long& destIdx)
  {
    if (originIdx >= m_nodes.size() || destIdx >= m_nodes.size())
    {
      ERROR_MSG("Invalid node indices " << originIdx << " and " << destIdx);
      abort();
    }
    if (m_nodes[originIdx] == 0 || m_nodes[destIdx] == 0)
    {
      ERROR_MSG("One of " << originIdx << " and " << destIdx << " does not exist as a node");
      abort();
    }
    if (originIdx == destIdx)
    {
      ERROR_MSG("Origin and destination nodes are the same!!");
      abort();
    }
    Edge* edgePtr = getEdge(edgeIdx);
    m_inEdges[edgePtr->m_toNodeIndex].first[edgePtr->m_toIndex] = -1;
    m_inEdges[edgePtr->m_toNodeIndex].second.push_front(edgePtr->m_toIndex);

    m_outEdges[edgePtr->m_fromNodeIndex].first[edgePtr->m_fromIndex] = -1;
    m_outEdges[edgePtr->m_fromNodeIndex].second.push_front(edgePtr->m_fromIndex);

    unsigned long inIdx = this->m_getFreeInAdacIdx(destIdx);
    unsigned long outIdx = this->m_getFreeOutAdacIdx(originIdx);

    edgePtr->m_fromIndex = outIdx;
    edgePtr->m_toIndex = inIdx;
    edgePtr->m_fromNodeIndex = originIdx;
    edgePtr->m_toNodeIndex = destIdx;

    m_outEdges[originIdx].first[outIdx] = edgeIdx;
    m_inEdges[destIdx].first[inIdx] = edgeIdx;

    return edgePtr;
  }

  void BaseGraph::removeNode(unsigned long nodeIdx)
  {
    if (nodeIdx >= m_nodes.size())
    {
      ERROR_MSG("Invalid node index " << nodeIdx);
      abort();
    }
    Node* node = m_nodes[nodeIdx];
    if (node == 0)
    {
      ERROR_MSG("Node at " << nodeIdx << " does not exist");
      abort();
    }
    delete node;
    m_nodes[nodeIdx] = 0;
    m_freeNodes.push_front(nodeIdx);

    pair<IndexVector, FreeList>& outEdges = m_outEdges[nodeIdx];
    for (unsigned long i = 0; i < outEdges.first.size(); i++)
    {
      const long& edgeIdx = outEdges.first[i];
      if (edgeIdx < 0)
      {
        continue;
      }
      Edge* edgePtr = m_edges[edgeIdx];
      m_inEdges[edgePtr->m_toNodeIndex].first[edgePtr->m_toIndex] = -1;
      m_inEdges[edgePtr->m_toNodeIndex].second.push_front(edgePtr->m_toIndex);

      m_edges[edgeIdx] = 0;
      m_freeEdges.push_front(edgeIdx);
      delete edgePtr;
    }
    m_outEdges[nodeIdx].first.resize(0);
    m_outEdges[nodeIdx].second.clear();

    pair<IndexVector, FreeList>& inEdges = m_inEdges[nodeIdx];
    for (unsigned long i = 0; i < inEdges.first.size(); i++)
    {
      const long& edgeIdx = inEdges.first[i];
      if (edgeIdx < 0)
      {
        continue;
      }
      Edge* edgePtr = m_edges[edgeIdx];

      m_outEdges[edgePtr->m_fromNodeIndex].first[edgePtr->m_fromIndex] = -1;
      m_outEdges[edgePtr->m_fromNodeIndex].second.push_front(edgePtr->m_fromIndex);

      m_edges[edgeIdx] = 0;
      m_freeEdges.push_front(edgeIdx);
      delete edgePtr;
    }
    m_inEdges[nodeIdx].first.resize(0);
    m_inEdges[nodeIdx].second.clear();
  }

  void BaseGraph::removeEdge(unsigned long edgeIdx)
  {
    if (edgeIdx >= m_edges.size())
    {
      ERROR_MSG("Invalid edge index " << edgeIdx);
      abort();
    }
    Edge* edgePtr = m_edges[edgeIdx];
    if (edgePtr == 0)
    {
      ERROR_MSG("Edge at " << edgeIdx << " does not exist");
      abort();
    }

    m_inEdges[edgePtr->m_toNodeIndex].first[edgePtr->m_toIndex] = -1;
    m_inEdges[edgePtr->m_toNodeIndex].second.push_front(edgePtr->m_toIndex);

    m_outEdges[edgePtr->m_fromNodeIndex].first[edgePtr->m_fromIndex] = -1;
    m_outEdges[edgePtr->m_fromNodeIndex].second.push_front(edgePtr->m_fromIndex);

    m_edges[edgeIdx] = 0;
    m_freeEdges.push_front(edgeIdx);
    delete edgePtr;
  }

  void BaseGraph::compress(vector<long>* outputNewNodeIdxs,
                           vector<long>* outputNewEdgeIdxs)
  {
    vector<long>* newNodeIdxs = (outputNewNodeIdxs == 0) ? new vector<long> ()
        : outputNewNodeIdxs;
    vector<long>* newEdgeIdxs = (outputNewEdgeIdxs == 0) ? new vector<long> ()
        : outputNewEdgeIdxs;

    newNodeIdxs->resize(m_nodes.size());
    newNodeIdxs->assign(m_nodes.size(), -1);
    newEdgeIdxs->resize(m_edges.size());
    newEdgeIdxs->assign(m_edges.size(), -1);

    unsigned long idxUse = 0;
    for (unsigned long i = 0; i < m_nodes.size(); i++)
    {
      if (m_nodes[i] == 0)
        continue;

      (*newNodeIdxs)[i] = idxUse;
      //DEBUG_MSG("Node " << i << " -> " << idxUse);
      m_nodes[idxUse] = m_nodes[i];
      m_nodes[idxUse]->m_index = idxUse;
      m_outEdges[idxUse].first = m_outEdges[i].first;
      m_outEdges[idxUse].second.clear();
      m_inEdges[idxUse].first = m_inEdges[i].first;
      m_inEdges[idxUse].second.clear();
      ++idxUse;
    }
    m_freeNodes.clear();
    m_nodes.resize(idxUse);

    idxUse = 0;
    for (unsigned long i = 0; i < m_edges.size(); i++)
    {
      if (m_edges[i] == 0)
        continue;

      (*newEdgeIdxs)[i] = idxUse;
      //DEBUG_MSG("Edge " << i << " -> " << idxUse);
      m_edges[idxUse] = m_edges[i];
      m_edges[idxUse]->m_fromNodeIndex
          = (*newNodeIdxs)[m_edges[idxUse]->m_fromNodeIndex];
      m_edges[idxUse]->m_toNodeIndex
          = (*newNodeIdxs)[m_edges[idxUse]->m_toNodeIndex];
      m_edges[idxUse]->m_index = idxUse;
      ++idxUse;
    }
    m_freeEdges.clear();
    m_edges.resize(idxUse);

    for (unsigned long i = 0; i < m_nodes.size(); i++)
    {
      idxUse = 0;
      for (unsigned long j = 0; j < m_outEdges[i].first.size(); j++)
      {
        if (m_outEdges[i].first[j] < 0)
          continue;

        long newEdgeIdx = (*newEdgeIdxs)[m_outEdges[i].first[j]];
        m_outEdges[i].first[idxUse] = newEdgeIdx;
        m_edges[newEdgeIdx]->m_fromIndex = idxUse;
        ++idxUse;
      }
      m_outEdges[i].first.resize(idxUse);

      idxUse = 0;
      for (unsigned long j = 0; j < m_inEdges[i].first.size(); j++)
      {
        if (m_inEdges[i].first[j] < 0)
          continue;

        long newEdgeIdx = (*newEdgeIdxs)[m_inEdges[i].first[j]];
        m_inEdges[i].first[idxUse] = newEdgeIdx;
        m_edges[newEdgeIdx]->m_toIndex = idxUse;
        ++idxUse;
      }
      m_inEdges[i].first.resize(idxUse);
    }

    m_outEdges.resize(m_nodes.size());
    m_inEdges.resize(m_nodes.size());

    if (outputNewNodeIdxs == 0)
      delete newNodeIdxs;

    if (outputNewEdgeIdxs == 0)
      delete newEdgeIdxs;
  }

  typedef sps::tuple<bool, double, Edge*> NodeTuple;

  class NodeDistCompare
  {
    bool reverse;
  public:
    NodeDistCompare(const bool& revparam = false)
    {
      reverse = revparam;
    }
    bool operator()(const NodeTuple& lhs, const NodeTuple& rhs) const
    {
      if (reverse)
        return (lhs.m1 >= rhs.m1);
      else
        return (lhs.m1 <= rhs.m1);
    }
  };

  typedef HashHeap<Node*, NodeTuple, NodeDistCompare> NodeDist;

  void BaseGraph::getLightestPaths(Node* source, Tree& outputPaths)
  {
    outputPaths.clear();
    if (!containsNode(source))
    {
      ERROR_MSG("Cannot find source \'" << source->toString() << "\'");
      return;
    }

    NodeDist distances(numNodes(), NodeDistCompare(false));

    const double maxDVal = numeric_limits<double>::max();

    NodeTuple defaultVal;
    defaultVal.m0 = false;
    defaultVal.m1 = maxDVal;
    defaultVal.m2 = 0;

    for (NodeSet::iterator nodeIt = m_nodes.begin(); nodeIt != m_nodes.end(); nodeIt++)
    {
      if (*nodeIt != 0)
      {
        distances.insert(*nodeIt, defaultVal);
      }
    }
    defaultVal.m0 = false;
    defaultVal.m1 = 0;
    defaultVal.m2 = 0;
    distances.update(source, defaultVal);

    while (distances.size() > 0)
    {
      pair<Node*, NodeTuple> heav = distances.popPair();

      Node* curNode = heav.first;
      double pathWeight = heav.second.m1;

      if (pathWeight == maxDVal)
      {
        break;
      }

      outputPaths[curNode] = pair<double, Edge*> (pathWeight, heav.second.m2);

      const IndexVector& outEdges = this->getOutEdges(curNode);

      for (unsigned long i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
        {
          continue;
        }
        Edge* edgePtr = m_edges[outEdges[i]];
        Node* toNode = m_nodes[edgePtr->getToNodeIndex()];

        if (distances.count(toNode) == 0)
        {
          continue;
        }

        NodeTuple distVis = distances[toNode];
        double curWeight = pathWeight + edgePtr->getWeight();

        if (curWeight < distVis.m1)
        {
          distVis.m1 = curWeight;
          distVis.m2 = edgePtr;
          distances.update(toNode, distVis);
        }
      }
    }
  }

  pair<bool, double> BaseGraph::getHeaviestPathDAG(Node* source,
                                                   Node* sink,
                                                   Path& outputPath)
  {
    if (!containsNode(source))
    {
      ERROR_MSG("Cannot find source \'" << source->toString() << "\'");
      return pair<bool, double> (false, 0);
    }
    if (!containsNode(sink))
    {
      ERROR_MSG("Cannot find sink \'" << sink->toString() << "\'");
      return pair<bool, double> (false, 0);
    }

    if (source->getIndex() == sink->getIndex())
    {
      WARN_MSG("Source and sink nodes are the same - \'" << source->toString() << "\'");
      return pair<bool, double> (false, 0);
    }
    list<Node*> topologOrder;
    if (!getTopologicalOrderingDAG(topologOrder))
    {
      ERROR_MSG("Graph contains a cycle, heaviest path may not work optimally");
    }

    vector<pair<double, Edge*> > longestDist(m_nodes.size());
    pair<double, Edge*> nextDist;
    const double minDVal = 0.0 - numeric_limits<double>::max();
    nextDist.first = minDVal;
    nextDist.second = 0;

    for (list<Node*>::const_iterator nIt = topologOrder.begin(); nIt
        != topologOrder.end(); nIt++)
    {
      longestDist[(*nIt)->getIndex()] = nextDist;
    }

    nextDist.first = 0;
    longestDist[source->getIndex()] = nextDist;

    for (list<Node*>::const_iterator nIt = topologOrder.begin(); nIt
        != topologOrder.end(); nIt++)
    {
      pair<double, Edge*>& curDist = longestDist[(*nIt)->getIndex()];
      if (curDist.first == minDVal)
      {
        continue;
      }

      const IndexVector& outEdges = this->getOutEdges(*nIt);
      for (unsigned long i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
          continue;

        Edge* edgePtr = m_edges[outEdges[i]];
        Node* toNode = m_nodes[edgePtr->getToNodeIndex()];
        double curWeight = curDist.first + edgePtr->getWeight();

        pair<double, Edge*>& neijDist = longestDist[toNode->getIndex()];

        if (curWeight > neijDist.first)
        {
          neijDist.first = curWeight;
          neijDist.second = edgePtr;
        }
      }
    }

    Node* nextNode = sink;
    list<Edge*> revPath;
    double totWeight = 0;

    do
    {
      pair<double, Edge*>& neijDist = longestDist[nextNode->getIndex()];
      revPath.push_back(neijDist.second);
      totWeight += neijDist.second->getWeight();
      nextNode = m_nodes[neijDist.second->getFromNodeIndex()];
    } while (nextNode->getIndex() != source->getIndex());

    outputPath.resize(revPath.size());
    unsigned int idxUse = 0;
    for (list<Edge*>::reverse_iterator eIt = revPath.rbegin(); eIt
        != revPath.rend(); eIt++)
    {
      outputPath[idxUse++] = *eIt;
    }

    return pair<bool, double> (true, totWeight);
  }

  bool BaseGraph::getTopologicalOrderingDAG(list<Node*>& outputOrder)
  {
    outputOrder.clear();
    list<Node*> sourceNodes;
    vector<unsigned int> incomingCount(m_nodes.size());
    unsigned int numEdges = 0;
    for (NodeSet::const_iterator nIt = m_nodes.begin(); nIt != m_nodes.end(); nIt++)
    {
      if (*nIt == 0)
        continue;

      unsigned int degree = getNumInEdges(*nIt);
      incomingCount[(*nIt)->getIndex()] = degree;
      numEdges += degree;
      if (degree == 0)
      {
        sourceNodes.push_back(*nIt);
      }
    }

    while (sourceNodes.size() > 0)
    {
      Node* curNode = sourceNodes.back();
      sourceNodes.pop_back();
      outputOrder.push_back(curNode);

      const IndexVector& outEdges = this->getOutEdges(curNode);
      for (unsigned long i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
          continue;

        Edge* edgePtr = m_edges[outEdges[i]];
        Node* toNode = m_nodes[edgePtr->getToNodeIndex()];

        incomingCount[toNode->getIndex()]--;
        numEdges--;
        if (incomingCount[toNode->getIndex()] == 0)
        {
          sourceNodes.push_back(toNode);
        }
      }
    }
    if (numEdges > 0)
    {
      WARN_MSG("Graph \'" << this->toString() << "\' has a cycle!");
      return false;
    }
    return true;
  }

  Node* BaseGraph::beginBFS(Node* source)
  {
    if (!containsNode(source))
    {
      return 0;
    }
    m_okTraverse = true;
    m_QueueIsStack = false;
    m_beginTraverse(source);
    return source;
  }
  
  Node* BaseGraph::beginBFS_Undirected(Node* source)
  {
    if (!containsNode(source))
    {
      return 0;
    }
    m_okTraverse = true;
    m_QueueIsStack = false;
    m_beginTraverse(source, true);
    return source;
  }

  Node* BaseGraph::beginDFS(Node* source)
  {
    if (!containsNode(source))
    {
      return 0;
    }
    m_okTraverse = true;
    m_QueueIsStack = true;
    m_beginTraverse(source);
    return source;
  }
  
  Node* BaseGraph::beginDFS_Undirected(Node* source)
  {
    if (!containsNode(source))
    {
      return 0;
    }
    m_okTraverse = true;
    m_QueueIsStack = true;
    m_beginTraverse(source, true);
    return source;
  }
  
  Node* BaseGraph::continueBFS_Undirected(Node* new_source)
  {
    if (!containsNode(new_source))
    {
      return 0;
    }
    
    if(m_visited.count(new_source) != 0){
        return 0;
    }
    
    m_okTraverse = true;
    m_QueueIsStack = false;
    m_continueTraverse(new_source, true);
    return new_source;
  }
  
  Node* BaseGraph::continueDFS_Undirected(Node* new_source)
  {
    if (!containsNode(new_source))
    {
      return 0;
    }
    
    if(m_visited.count(new_source) != 0){
        return 0;
    }
    
    m_okTraverse = true;
    m_QueueIsStack = true;
    m_continueTraverse(new_source, true);
    return new_source;
  }

  Node* BaseGraph::nextBFSDFS(void)
  {
    if (m_queueStack.size() == 0 || !m_okTraverse)
    {
      m_visited.clear();
      return 0;
    }
    Edge* next;
    if (m_QueueIsStack)
    {
      next = m_queueStack.back();
      m_queueStack.pop_back();
    }
    else
    {
      next = m_queueStack.front();
      m_queueStack.pop_front();
    }
    Node* to = m_nodes[next->getToNodeIndex()];
    const IndexVector& outEdges = getOutEdges(to);

    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      Edge* edgePtr = m_edges[outEdges[i]];
      Node* toNode = m_nodes[edgePtr->getToNodeIndex()];
      if (m_visited.count(toNode) == 0)
      {
        m_queueStack.push_back(edgePtr);
        m_visited.insert(toNode);
      }
    }
    return to;
  }
  
  
  Node* BaseGraph::nextBFSDFS_Undirected(Edge * &traversed_edge)
  {
    if (m_queueStack.size() == 0 || !m_okTraverse)
    {
      //m_visited.clear();
      return 0;
    }
    Edge* next;
    int traverse_direction;
    if (m_QueueIsStack)
    {
      next = m_queueStack.back();
      m_queueStack.pop_back();
      traverse_direction = m_queueStack_undirected_traversedirection.back();
      m_queueStack_undirected_traversedirection.pop_back();
    }
    else
    {
      next = m_queueStack.front();
      m_queueStack.pop_front();
      traverse_direction = m_queueStack_undirected_traversedirection.front();
      m_queueStack_undirected_traversedirection.pop_front();
    }
    Node* to;
    traversed_edge = next;
    if(traverse_direction == 0)
        to = m_nodes[next->getToNodeIndex()];
    else
        to = m_nodes[next->getFromNodeIndex()];
    const IndexVector& outEdges = getOutEdges(to);
    const IndexVector& inEdges = getInEdges(to);

    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      Edge* edgePtr = m_edges[outEdges[i]];
      Node* toNode = m_nodes[edgePtr->getToNodeIndex()];
      if (m_visited.count(toNode) == 0)
      {
        m_queueStack.push_back(edgePtr);
        m_visited.insert(toNode);
        m_queueStack_undirected_traversedirection.push_back(0);
      }
    }
    
    for (unsigned long i = 0; i < inEdges.size(); i++)
    {
      if (inEdges[i] < 0)
        continue;

      Edge* edgePtr = m_edges[inEdges[i]];
      Node* fromNode = m_nodes[edgePtr->getFromNodeIndex()];
      if (m_visited.count(fromNode) == 0)
      {
        m_queueStack.push_back(edgePtr);
        m_visited.insert(fromNode);
        m_queueStack_undirected_traversedirection.push_back(1);
      }
    }
    
    return to;
  }

  void BaseGraph::m_beginTraverse(Node* source, bool undirected)
  {
    m_queueStack.clear();
    m_queueStack_undirected_traversedirection.clear();
    m_visited.clear();
    

    m_continueTraverse(source, undirected);
  }
  
  void BaseGraph::m_continueTraverse(Node* source, bool undirected)
  {
    m_visited.insert(source);
    
    const IndexVector& outEdges = getOutEdges(source);
    const IndexVector& inEdges = getInEdges(source);

    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      Edge* edgePtr = m_edges[outEdges[i]];
      Node* toNode = m_nodes[edgePtr->getToNodeIndex()];
      if (m_visited.count(toNode) == 0)
      {
        m_queueStack.push_back(edgePtr);
        m_visited.insert(toNode);
        m_queueStack_undirected_traversedirection.push_back(0);
      }
    }
    
    if(undirected){
        for (unsigned long i = 0; i < inEdges.size(); i++)
        {
            if (inEdges[i] < 0)
                continue;

            Edge* edgePtr = m_edges[inEdges[i]];
            Node* fromNode = m_nodes[edgePtr->getFromNodeIndex()];
            if (m_visited.count(fromNode) == 0)
            {
                m_queueStack.push_back(edgePtr);
                m_visited.insert(fromNode);
                m_queueStack_undirected_traversedirection.push_back(1);
            }
        }
    }
  }
  
  void BaseGraph::clear_traversal_datastructures(){
        m_queueStack.clear();
        m_queueStack_undirected_traversedirection.clear();
        m_visited.clear();
    }       
  
  void BaseGraph::m_clearNodes(void)
  {
    for (unsigned long i = 0; i < m_edges.size(); i++)
    {
      if (m_edges[i] != 0)
      {
        delete m_edges[i];
      }
    }
    for (unsigned long i = 0; i < m_nodes.size(); i++)
    {
      if (m_nodes[i] != 0)
      {
        delete m_nodes[i];
        m_nodes[i] = 0;
        m_outEdges[i].first.resize(0);
        m_inEdges[i].first.resize(0);
        m_outEdges[i].second.clear();
        m_inEdges[i].second.clear();
      }
    }
  }

  unsigned long BaseGraph::m_getFreeNodeIdx(void)
  {
    unsigned long newIdx;
    if (m_freeNodes.size() == 0)
    {
      newIdx = m_nodes.size();
      m_nodes.resize(newIdx + 1);
      m_inEdges.resize(newIdx + 1);
      m_outEdges.resize(newIdx + 1);
    }
    else
    {
      newIdx = m_freeNodes.front();
      m_freeNodes.pop_front();
    }
    return newIdx;
  }

  unsigned long BaseGraph::m_getFreeEdgeIdx(void)
  {
    unsigned long newIdx;
    if (m_freeEdges.size() == 0)
    {
      newIdx = m_edges.size();
      m_edges.resize(newIdx + 1);
    }
    else
    {
      newIdx = m_freeEdges.front();
      m_freeEdges.pop_front();
    }
    return newIdx;
  }

  unsigned long BaseGraph::m_getFreeInAdacIdx(unsigned long nodeIdx)
  {
    unsigned long newIdx;
    if (m_inEdges[nodeIdx].second.size() == 0)
    {
      newIdx = m_inEdges[nodeIdx].first.size();
      m_inEdges[nodeIdx].first.resize(newIdx + 1);
    }
    else
    {
      newIdx = m_inEdges[nodeIdx].second.front();
      m_inEdges[nodeIdx].second.pop_front();
    }
    return newIdx;
  }

  unsigned long BaseGraph::m_getFreeOutAdacIdx(unsigned long nodeIdx)
  {
    unsigned long newIdx;
    if (m_outEdges[nodeIdx].second.size() == 0)
    {
      newIdx = m_outEdges[nodeIdx].first.size();
      m_outEdges[nodeIdx].first.resize(newIdx + 1);
    }
    else
    {
      newIdx = m_outEdges[nodeIdx].second.front();
      m_outEdges[nodeIdx].second.pop_front();
    }
    return newIdx;
  }
}
