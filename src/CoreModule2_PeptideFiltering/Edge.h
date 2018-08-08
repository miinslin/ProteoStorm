/*
 * Edge.h
 *
 *  Created on: Feb 7, 2012
 *      Author: aguthals
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "Node.h"
#include "Logger.h"

using namespace std;

namespace specnets
{

  class Node;

  extern class EdgeMismatchException : public exception
  {
  public:
    virtual const char* what() const throw ()
    {
      return "The edge's source/sink node(s) do not match their expected values";
    }
  } EMExcep;

  extern class EdgeNotFoundException : public exception
  {
  public:
    virtual const char* what() const throw ()
    {
      return "Edge not found";
    }
  } ENFExcep;

  class Edge
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    // Default constructor. Should always be called by derived class constructors
    Edge(void);

    // Copy constructor. Should always be called by derived class constructors
    Edge(const Edge& other);

    virtual ~Edge();

    virtual string toString(void) const;

    virtual string getGraphvizLabel(void) const;

    Edge &operator=(const Edge &other);

    // copies another edge
    virtual void copy(Edge& other);

    // Sets the edge weight
    inline virtual void setWeight(const double& newWeight)
    {
      weight = newWeight;
    }

    // Adds to the edge weight
    inline virtual void addWeight(const double& addWeight)
    {
      weight += addWeight;
    }

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

    inline const double& getWeight(void) const
    {
      return weight;
    }

    /**
     * Less than operator.
     */
    inline bool operator<(const Edge& other) const
    {
      return weight < other.getWeight();
    }

    /**
     * Greater than operator.
     */
    inline bool operator>(const Edge& other) const
    {
      return weight > other.getWeight();
    }

    /**
     * Equals operator.
     */
    inline bool operator==(const Edge& other) const
    {
      return weight == other.getWeight();
    }

    /**
     * Not equals operator.
     */
    inline bool operator!=(const Edge& other) const
    {
      return weight != other.getWeight();
    }

    /**
     * <= operator
     */
    inline bool operator<=(const Edge& other) const
    {
      return weight <= other.getWeight();
    }

    /**
     * >= operator
     */
    inline bool operator>=(const Edge& other) const
    {
      return weight >= other.getWeight();
    }

    inline const long& getIndex(void) const
    {
      return m_index;
    }
    inline const long& getFromIndex(void) const
    {
      return m_fromIndex;
    }
    inline const long& getToIndex(void) const
    {
      return m_toIndex;
    }
    inline const long& getFromNodeIndex(void) const
    {
      return m_fromNodeIndex;
    }
    inline const long& getToNodeIndex(void) const
    {
      return m_toNodeIndex;
    }

  protected:
    // Edge weight
    double weight;

  private:
    // Source node
    long m_fromNodeIndex;

    // Sink node
    long m_toNodeIndex;

    // This edge's index in BaseGraph
    long m_index;

    // This edge's index in the from Node's adjacency vector in BaseGraph
    long m_fromIndex;

    // This edge's index in the to Node's adjacency vector in BaseGraph
    long m_toIndex;

    friend class BaseGraph;

  };

  class EdgePtrCompare
  {
    bool m_reverse;
  public:
    EdgePtrCompare(const bool& revparam = false)
    {
      m_reverse = revparam;
    }
    bool operator()(const Edge* lhs, const Edge* rhs) const
    {
      if (m_reverse)
        return (lhs->getWeight() >= rhs->getWeight());
      else
        return (lhs->getWeight() <= rhs->getWeight());
    }
  };

}

#endif /* EDGE_H_ */
