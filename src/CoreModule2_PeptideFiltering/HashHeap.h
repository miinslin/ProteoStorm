/*
 * HashHeap.h
 *
 *  Created on: Feb 12, 2012
 *      Author: aguthals
 *
 *  Description: HashHeaps are associative containers that store elements formed by key-value
 *    pairs and keep a weak-ordering of the mapped values. This allows for fast retrieval of individual
 *    elements based on their keys as well as log(n) retrieval of elements with the maximum/minimum value. Its
 *    purpose was motivated by the inability of STL's priority_queue to efficiently increase/decrease values in
 *    place. Thus, this is essentially a priority queue that supports efficient modification of heap values. The
 *    purpose of the hash table (tr1::unordered_map) is to map elements' keys to their position in the Heap.
 */

#ifndef HASHHEAP_H_
#define HASHHEAP_H_

#include <vector>
#include <math.h>

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

#include "Logger.h"
#include "utils.h"

using namespace std;

namespace specnets
{

  /**
   * Exception thrown when an input key is not found in the hash table
   */
  class KeyNotFoundException : public exception
  {
    virtual const char* what() const throw ()
    {
      return "Key not found in hash table";
    }
  } NoKeyExcep;

  /**
   * Key - Type of the key values. Each element is uniquely identified by its key value via a hash table.
   * T - Type of the ordered mapped values. Elements are ordered by their mapped values in the Heap.
   * Compare - Comparison class: A class such that the expression comp(a,b), where comp is an object of
   *   this class and a and b are elements of the container, returns true if a is to be placed earlier
   *   than b in a strict weak ordering operation.
   * Hash - A unary function object type that takes an object of type key type as argument and returns a
   *   unique value of type size_t based on it. See http://www.cplusplus.com/reference/stl/unordered_map/
   * Pred - A binary predicate that takes two arguments of the key type and returns a bool.
   *   See http://www.cplusplus.com/reference/stl/unordered_map/
   */
  template<class Key, class T = Key, class Compare = less<T> ,
      class Hash = tr1::hash<Key>, class Pred = equal_to<Key> > class HashHeap
  {

  protected:

    typedef vector<pair<Key, T> > Array;
    typedef tr1::unordered_map<Key, size_t, Hash, Pred> IndexMap;

    // Heap stored in a vector as a balanced binary tree
    Array* Tree;

    // Hash table for referencing elements' position in the Heap
    IndexMap* IdxMap;

    // Heap index of the next element to be added. This is also the number of elements in the Heap
    size_t nextIdx;

    // Comparator object for ordering mapped values
    const Compare& Comp;

    // tr1::unordered_map template parameters
    const Hash& HashObj;
    const Pred& PredObj;

  public:
    /**
     * Initializes heap size to 0.
     *@param x instance of Compare to use when ordering values
     *@param h instance of Hash fo hashing keys to their values
     *@param r instance of Pred for determining equality of keys
     */
    HashHeap(const Compare& x = Compare(),
             const Hash& h = Hash(),
             const Pred& p = Pred()) :
      Tree(0x0), Comp(x), IdxMap(0x0), nextIdx(0), HashObj(h), PredObj(p)
    {
      initialize(0);
    }

    /**
     * Initializes heap size.
     *@param sz Number of elements to make room for - NOTE: This does not insert any elements into the HashHeap
     *@param x instance of Compare to use when ordering values
     *@param h instance of Hash fo hashing keys to their values
     *@param r instance of Pred for determining equality of keys
     */
    HashHeap(const size_t& sz, const Compare& x = Compare(), const Hash& h =
        Hash(), const Pred& p = Pred()) :
      Tree(0x0), Comp(x), IdxMap(0x0), nextIdx(0), HashObj(h), PredObj(p)
    {
      initialize(sz);
    }

    /**
     * Deallocates memory for Tree and IdxMap
     */
    ~HashHeap()
    {
      if (Tree)
      {
        delete Tree;
      }
      if (IdxMap)
      {
        delete IdxMap;
      }
    }

    /**
     * Initializes data structures. If the HashHeap was previously initialized, this will remove all the elements
     *@param sz Number of elements to make room for - NOTE: This does not insert any elements into the HashHeap
     */
    void initialize(const size_t& sz)
    {
      if (Tree)
      {
        Tree->resize(sz);
      }
      else
      {
        Tree = new Array(sz);
      }

      if (IdxMap)
      {
        IdxMap->clear();
        IdxMap->rehash(sz);
      }
      else
      {
        IdxMap = new IndexMap(sz, HashObj, PredObj);
      }

      nextIdx = 0;
    }

    /**
     * Returns the number of elements stored on the heap and referenced by the hash table
     */
    size_t size(void) const
    {
      return nextIdx;
    }

    /**
     * Returns the maximum number of elements that the HashHeap has room for.
     */
    size_t capacity(void) const
    {
      return Tree->size();
    }

    /**
     * Decreases capacity of the HashHeap to its size.
     */
    void compress(void)
    {
      Tree->resize(nextIdx);
      IdxMap->rehash(nextIdx);
    }

    /**
     * Increases the capacity of the HashHeap
     *@param deltaSz positive change in capacity
     */
    void expand(const size_t& deltaSz)
    {
      const size_t newSize = Tree->size() + deltaSz;
      Tree->resize(newSize);
      IdxMap->rehash(newSize);
    }

    /**
     * Returns the number of elements matching the key. Constant average running time.
     */
    size_t count(const Key& key) const
    {
      return IdxMap->count(key);
    }

    /**
     * Returns a reference to the value mapped by the key. Constant average running time.
     */
    T & operator[](const Key& key) const
    {
      if (IdxMap->count(key) == 0)
      {
        ERROR_MSG("Key \'" << key << "\' not found");
        throw NoKeyExcep;
      }
      return (*Tree)[(*IdxMap)[key]].second;
    }

    /**
     * Inserts a new element into the HashHeap. If an element with a matching key is already stored, its value is updated with the input value.
     *   If the key is unique and the HashHeap is already filled to its capacity, the capacity is increased by one to make room.
     *   Logarithmic running time when the new element is the new maximum, constant running time when the new element is the minimum value
     *@param key Key instance to reference new element by
     *@param item T instance to order new element by
     *@return
     */
    void insert(const Key& key, const T& item)
    {
      if (IdxMap->count(key) > 0)
      {
        update(key, item);
        return;
      }

      const size_t newIdx = nextIdx++;

      if (newIdx >= IdxMap->size())
      {
        IdxMap->rehash(newIdx + 1);
        Tree->resize(newIdx + 1);
      }

      (*IdxMap)[key] = newIdx;
      (*Tree)[newIdx].first = key;
      (*Tree)[newIdx].second = item;

      m_increaseKey(newIdx);
    }

    /**
     * Returns a reference to the maximum mapped value in the HashHeap. Constant running time.
     */
    T& top(void) const
    {
      if (size() == 0)
      {
        throw NoKeyExcep;
      }
      return (*Tree)[0].second;

    }

    /**
     * Returns a reference to the maximum element in the HashHeap (as ordered by mapped values). Constant running time.
     */
    pair<Key, T>& topPair(void) const
    {
      if (size() == 0)
      {
        throw NoKeyExcep;
      }
      return (*Tree)[0];
    }

    /**
     * Removes and returns the maximum mapped value in the HashHeap. Logarithmic running time in the worst case.
     */
    T pop(void)
    {
      if (size() == 0)
      {
        throw NoKeyExcep;
      }

      pair<Key, T> root = popPair();

      return T(root.second);
    }

    /**
     * Removes and returns the maximum element in the HashHeap. Logarithmic running time in the worst case.
     */
    pair<Key, T> popPair(void)
    {
      if (size() == 0)
      {
        throw NoKeyExcep;
      }

      pair<Key, T> root((*Tree)[0].first, (*Tree)[0].second);

      if (size() == 1)
      {
        trim();
        return pair<Key, T> (root.first, root.second);
      }

      swap(nextIdx - 1, 0);
      trim();
      m_decreaseKey(0);
      return pair<Key, T> (root.first, root.second);
    }

    /**
     * Removes an element identified by key from the HashHeap. Logarithmic running time in the worst case.
     *@param key Key of element to remove
     *@return
     */
    void erase(const Key& key)
    {
      if (IdxMap->count(key) == 0)
      {
        throw NoKeyExcep;
      }
      if (size() == 1)
      {
        trim();
        return;
      }
      const size_t keyIdx = (*IdxMap)[key];
      const T& oldVal = (*Tree)[keyIdx].second;
      const T& newVal = (*Tree)[nextIdx - 1].second;
      bool inc = Comp(newVal, oldVal);
      swap(nextIdx - 1, keyIdx);
      trim();
      if (inc)
      {
        m_increaseKey(keyIdx);
      }
      else
      {
        m_decreaseKey(keyIdx);
      }
    }

    /**
     * If the value mapped by key has increased, this will update the HashHeap accordingly. Logarithmic worst case running time.
     *@param key Key of element whose mapped value has increased
     *@return
     */
    void increase(const Key& key)
    {

      if (IdxMap->count(key) == 0)
      {
        throw NoKeyExcep;
      }
      m_increaseKey((*IdxMap)[key]);
    }

    /**
     * If the value mapped by key has decreased, this will update the HashHeap accordingly. Logarithmic worst case running time.
     *@param key Key of element whose mapped value has decreased
     *@return
     */
    void decrease(const Key& key)
    {
      if (IdxMap->count(key) == 0)
      {
        throw NoKeyExcep;
      }
      m_decreaseKey((*IdxMap)[key]);
    }

    /**
     * Updates the mapped value of a element. This effectively calls increase() if the element's value has increased and
     *   decrease() if the element's value has decreased. Logarithmic worst case running time.
     *@param key
     *@param newValue
     *@return
     */
    void update(const Key& key, const T& newValue)
    {
      if (IdxMap->count(key) == 0)
      {
        throw NoKeyExcep;
      }
      const size_t& keyIdx = (*IdxMap)[key];
      const T oldVal = (*Tree)[keyIdx].second;
      (*Tree)[keyIdx].second = newValue;

      //DEBUG_VAR(Comp(oldVal, newValue));
      //DEBUG_VAR(Comp(newValue, oldVal));

      if (Comp(newValue, oldVal))
      {
        m_increaseKey(keyIdx);
      }
      else
      {
        m_decreaseKey(keyIdx);
      }
    }

  private:
    // internal decrease()
    void m_decreaseKey(const size_t& inIdx)
    {
      size_t atIdx = inIdx;

      size_t lChild = getLChildIdx(atIdx);
      size_t rChild = getRChildIdx(atIdx);
      size_t swapChild;

      bool updateLeft = lChild < nextIdx && (!Comp((*Tree)[atIdx].second,
                                                   (*Tree)[lChild].second));
      bool updateRight = rChild < nextIdx && (!Comp((*Tree)[atIdx].second,
                                                    (*Tree)[rChild].second));

      while (updateLeft || updateRight)
      {
        if (lChild >= nextIdx)
        {
          swapChild = rChild;
        }
        else if (rChild >= nextIdx)
        {
          swapChild = lChild;
        }
        else
        {
          swapChild = (Comp((*Tree)[rChild].second, (*Tree)[lChild].second))
              ? rChild : lChild;
        }
        swap(swapChild, atIdx);
        atIdx = swapChild;
        lChild = getLChildIdx(atIdx);
        rChild = getRChildIdx(atIdx);
        updateLeft = lChild < nextIdx && (!Comp((*Tree)[atIdx].second,
                                                (*Tree)[lChild].second));
        updateRight = rChild < nextIdx && (!Comp((*Tree)[atIdx].second,
                                                 (*Tree)[rChild].second));
      }
    }

    // internal increase()
    void m_increaseKey(const size_t& inIdx)
    {
      size_t atIdx = inIdx;
      size_t parIdx = getParIdx(atIdx);

      while (parIdx != atIdx && (!Comp((*Tree)[parIdx].second,
                                       (*Tree)[atIdx].second)))
      {
        swap(atIdx, parIdx);

        atIdx = parIdx;
        parIdx = getParIdx(atIdx);
      }
    }

    // removes the last element in the heap
    void trim(void)
    {
      const Key& lastKey = (*Tree)[nextIdx - 1].first;
      IdxMap->erase(lastKey);
      nextIdx--;
    }

    // swaps two heap nodes and update hash mappings
    void swap(const size_t& childIdx, const size_t& parIdx)
    {

      const Key childID = (*Tree)[childIdx].first;
      const Key parID = (*Tree)[parIdx].first;
      const T parItem = (*Tree)[parIdx].second;

      (*Tree)[parIdx].first = childID;
      (*Tree)[parIdx].second = (*Tree)[childIdx].second;

      (*Tree)[childIdx].first = parID;
      (*Tree)[childIdx].second = parItem;

      (*IdxMap)[parID] = childIdx;
      (*IdxMap)[childID] = parIdx;
    }

    // gets parent index of child in the tree vector
    static inline size_t getParIdx(const size_t& childIdx)
    {
      return (childIdx == 0) ? 0 : floatToInt(floor(((double)(childIdx - 1))
          / 2.0));
    }

    // gets left child index of parent in the tree vector
    static inline size_t getLChildIdx(const size_t& parIdx)
    {
      return (parIdx * 2) + 1;
    }

    // gets right child index of parent in the tree vector
    static inline size_t getRChildIdx(const size_t& parIdx)
    {
      return (parIdx * 2) + 2;
    }

  };
}

#endif /* HASHHEAP_H_ */
