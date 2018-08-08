#ifndef TWOVALUES_H
#define TWOVALUES_H

#include <iostream>

/**
* TODO: add description
*/
template<class T> class TwoValues {
public:
    /**
    * TODO: add description
    */

    TwoValues():
      values((T) 0, (T) 0)
    {
    }

    /**
    * TODO: add description
    */
    TwoValues(T a, T b) {
        values.first = a;
        values.second = b;
    }

    /**
    * TODO: add description
    *
    *@param other
    */
    TwoValues(const TwoValues<T> &other) {
        values.first = other.values.first;
        values.second = other.values.second;
    }

    /**
    * TODO: add description
    *
    *@param other
    *@return
    */
    TwoValues<T> &operator=(const TwoValues<T> &other) {
        values.first = other.values.first;
        values.second = other.values.second;
        return (*this);
    }

    /**
    * TODO: add description
    *
    *@param a
    *@param b
    *@return
    */
    TwoValues<T> &set(T a, T b) {
        values.first = a;
        values.second = b;
        return (*this);
    }

    /**
    * TODO: add description
    *
    *@param i
    *@return
    */
    T &operator[](int i) {
      switch(i)
      {
      case 0:
        return values.first;
      case 1:
        return values.second;
      default:
        std::cerr << "Invalid input value for i: " << i << std::endl;
        return values.first;
      }
    }

    const T &operator[](int i) const {
      switch(i)
      {
      case 0:
        return values.first;
      case 1:
        return values.second;
      default:
        std::cerr << "Invalid input value for i: " << i << std::endl;
        return values.first;
      }
    }

    /**
    * TODO: add description
    *
    *@param other
    *@return
    */
    bool operator<(const TwoValues<T> &other) {
        if (values.first < other.values.first or (values.first == other.values.first
                and values.second < other.values.second))
            return true;
        return false;
    }

    /**
    * TODO: add description
    *
    *@param other
    *@return
    */
    bool operator==(const TwoValues<T> &other) {
        return values.first == other.values.first and values.second == other.values.second;
    }
    //    void operator<<(ostream &out) { out<<"("<<values[0]<<","<<values[1]<<")"; }
public:
  std::pair<T,T > values;
};

/**
* TODO: add description
*/
template<class T> bool operator<(const TwoValues<T> &a, const TwoValues<T> &b) {
    if (a[0] < b[0] or (a[0] == b[0]
            and a[1] < b[1]))
        return true;
    return false;
}



#endif
