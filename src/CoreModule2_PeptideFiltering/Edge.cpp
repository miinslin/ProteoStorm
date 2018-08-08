/*
 * Edge.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: aguthals
 */

#include "Edge.h"

using namespace std;

namespace specnets
{

  const unsigned short Edge::BIN_VERSION = 1;
  const unsigned short Edge::BIN_SUBVERSION = 1;

  const string Edge::BIN_VERSION_ID = "Edge_binVersion";
  const string Edge::BIN_SUBVERSION_ID = "Edge_binSubVersion";

  Edge::Edge(void) :
    weight(0), m_index(-1), m_fromIndex(-1), m_toIndex(-1),
        m_fromNodeIndex(-1), m_toNodeIndex(-1)
  {
    //DEBUG_MSG("Constructor for edge <" << this << ">");
  }

  Edge::Edge(const Edge& other) :
    weight(other.weight), m_index(other.m_index),
        m_fromIndex(other.m_fromIndex), m_toIndex(other.m_toIndex),
        m_fromNodeIndex(other.m_fromNodeIndex),
        m_toNodeIndex(other.m_toNodeIndex)
  {
    //DEBUG_MSG("Copy constructor for " << this->toString());
  }

  Edge::~Edge()
  {
    //DEBUG_MSG("Destructor for " << this->toString());
  }

  string Edge::toString(void) const
  {
    ostringstream out;
    out << "Edge " << m_fromNodeIndex << " -> " << m_toNodeIndex << " at "
        << m_index;
    return out.str();
  }

  string Edge::getGraphvizLabel(void) const
  {
    ostringstream out;
    out << weight;
    return out.str();
  }

  Edge & Edge::operator=(const Edge &other)
  {
    if (this == &other)
    {
      return *this;
    }

    this->m_fromNodeIndex = other.m_fromNodeIndex;
    this->m_toNodeIndex = other.m_toNodeIndex;
    this->weight = other.getWeight();
    this->m_index = other.getIndex();
    this->m_fromIndex = other.m_fromIndex;
    this->m_toIndex = other.m_toIndex;

    return *this;
  }

  void Edge::copy(Edge& other)
  {
    this->operator =(other);
  }

  bool Edge::saveToBinaryStream(FILE* fp) const
  {
    unsigned int count;

    count = fwrite(&m_fromNodeIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the origin ID");
      return false;
    }
    count = fwrite(&m_toNodeIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the destination ID");
      return false;
    }
    count = fwrite(&m_index, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the index");
      return false;
    }
    count = fwrite(&m_fromIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the from index");
      return false;
    }
    count = fwrite(&m_toIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the to index");
      return false;
    }
    count = fwrite(&weight, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the weight");
      return false;
    }
    return true;
  }

  bool Edge::loadFromBinaryStream(FILE* fp,
                                  map<string, unsigned short>& versions)
  {

    unsigned int count;

    unsigned short edgeVersion = versions[BIN_VERSION_ID];
    unsigned short edgeSubVersion = versions[BIN_SUBVERSION_ID];

    count = fread(&m_fromNodeIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the origin ID");
      return false;
    }
    count = fread(&m_toNodeIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the destination ID");
      return false;
    }
    count = fread(&m_index, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the index");
      return false;
    }
    count = fread(&m_fromIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the from index");
      return false;
    }
    count = fread(&m_toIndex, sizeof(long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the to index");
      return false;
    }
    count = fread(&weight, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the weight");
      return false;
    }
    return true;
  }

}
