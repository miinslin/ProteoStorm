/*
 * Node.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: aguthals
 */

#include "Node.h"

using namespace std;

namespace specnets
{

  const unsigned short Node::BIN_VERSION = 1;
  const unsigned short Node::BIN_SUBVERSION = 1;

  const string Node::BIN_VERSION_ID = "Node_binVersion";
  const string Node::BIN_SUBVERSION_ID = "Node_binSubVersion";

  Node::Node() :
    m_index(-1)
  {
    this->initialize();
    //DEBUG_MSG("Constructor for node <" << this << ">");
  }

  Node::Node(const Node& other) :
    m_index(other.getIndex())
  {
    this->initialize();
    //DEBUG_MSG("Copy constructor for node <" << this << ">");
  }

  Node::~Node(void)
  {
    //DEBUG_MSG("Destructor for node " << getID() << " <" << this << ">");
  }

  string Node::toString() const
  {
    ostringstream out;
    out << m_index << " at <" << this << ">";
    return out.str();
  }

  string Node::getGraphvizLabel(void) const
  {
    ostringstream out;
    out << m_index;
    return out.str();
  }

  string Node::getGraphvizFillColor(void) const
  {
    return "";
  }

  Node & Node::operator=(const Node &other)
  {
    if (this == &other)
    {
      return *this;
    }
    m_index = other.getIndex();
    return *this;
  }

  void Node::copy(Node& other)
  {
    this->operator =(other);
  }

  bool Node::saveToBinaryStream(FILE* fp) const
  {
    unsigned int count;

    count = fwrite(&m_index, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the node ID");
      return false;
    }
    return true;
  }

  bool Node::loadFromBinaryStream(FILE* fp,
                                  map<string, unsigned short>& versions)
  {

    unsigned short nodeVersion = versions[BIN_VERSION_ID];
    unsigned short nodeSubVersion = versions[BIN_SUBVERSION_ID];

    unsigned int count = fread(&m_index, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the node ID");
      return false;
    }
    return true;
  }

  void Node::initialize(void)
  {
  }

}
