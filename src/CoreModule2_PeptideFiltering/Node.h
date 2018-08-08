/*
 * Node.h
 *
 *  Created on: Feb 7, 2012
 *      Author: aguthals
 */

#ifndef NODE_H_
#define NODE_H_

#include <map>
#include <list>

#include "Edge.h"
#include "Logger.h"

using namespace std;

namespace specnets
{
  class Edge;
  class Node;

  extern class NodeNotFoundException : public exception
  {
  public:
    virtual const char* what() const throw ()
    {
      return "Node not found";
    }
  } NNFExcep;

  class Node
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    Node();

    Node(const Node& other);

    virtual ~Node(void);

    virtual string toString(void) const;

    virtual string getGraphvizLabel(void) const;

    virtual string getGraphvizFillColor(void) const;

    Node &operator=(const Node &other);

    virtual void copy(Node& other);

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

    void initialize(void);

    inline const long& getIndex(void) const
    {
      return m_index;
    }

  private:

    // set by BaseGraph for indexing
    long m_index;

    friend class BaseGraph;
  };

}

#endif /* NODE_H_ */
