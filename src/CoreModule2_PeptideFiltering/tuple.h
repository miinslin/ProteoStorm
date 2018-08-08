/**
    @file tuple.h

    @note
    Copyright 2006, The Regents of the University of California
    All Rights Reserved

    Permission to use, copy, modify and distribute any part of this
    program for educational, research and non-profit purposes, without fee,
    and without a written agreement is hereby granted, provided that the
    above copyright notice, this paragraph and the following three paragraphs
    appear in all copies.

    Those desiring to incorporate this work into commercial
    products or use for commercial purposes should contact the Technology
    Transfer & Intellectual Property Services, University of California,
    San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
    Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.

    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
    FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
    INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
    IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
    OF SUCH DAMAGE.

    THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
    OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
    ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO
    REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
    EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
    MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
    THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
*/


#ifndef TUPLE_H
#define TUPLE_H


namespace sps
{

template <typename T0, typename T1 = void, typename T2 = void, typename T3 = void, typename T4 = void, typename T5 = void, typename T6 = void, typename T7 = void>
  struct tuple;

template <typename T0>
  struct tuple<T0, void, void, void, void, void, void, void>
  {
    T0 m0;
  };

template <typename T0, typename T1>
  struct tuple<T0, T1, void, void, void, void, void, void> : tuple<T0>
  {
    T1 m1;
  };

template <typename T0, typename T1, typename T2>
  struct tuple<T0, T1, T2, void, void, void, void, void> : tuple<T0, T1>
  {
    T2 m2;
  };

template <typename T0, typename T1, typename T2, typename T3>
  struct tuple<T0, T1, T2, T3, void, void, void, void> : tuple<T0, T1, T2>
  {
    T3 m3;
  };

template <typename T0, typename T1, typename T2, typename T3, typename T4>
  struct tuple<T0, T1, T2, T3, T4, void, void, void> : tuple<T0, T1, T2, T3>
  {
    T4 m4;
  };

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
  struct tuple<T0, T1, T2, T3, T4, T5, void, void> : tuple<T0, T1, T2, T3, T4>
  {
    T5 m5;
  };

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  struct tuple<T0, T1, T2, T3, T4, T5, T6, void> : tuple<T0, T1, T2, T3, T4, T5>
  {
    T6 m6;
  };

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  struct tuple : tuple<T0, T1, T2, T3, T4, T5, T6>
  {
    T7 m7;
  };

template <typename T0>
  tuple<T0> make_tuple(const T0 & m0)
  {
    tuple<T0> t;
    t.m0 = m0;
    return t;
  }

template <typename T0, typename T1>
  tuple<T0, T1> make_tuple(const T0 & m0, const T1 & m1)
  {
    tuple<T0, T1> t;
    t.m0 = m0;
    t.m1 = m1;
    return t;
  }

template <typename T0, typename T1, typename T2>
  tuple<T0, T1, T2> make_tuple(const T0 & m0, const T1 & m1, const T2 & m2)
  {
    tuple<T0, T1, T2> t;
    t.m0 = m0;
    t.m1 = m1;
    t.m2 = m2;
    return t;
  }

template <typename T0, typename T1, typename T2, typename T3>
  tuple<T0, T1, T2, T3> make_tuple(const T0 & m0, const T1 & m1, const T2 & m2, const T3 & m3)
  {
    tuple<T0, T1, T2, T3> t;
    t.m0 = m0;
    t.m1 = m1;
    t.m2 = m2;
    t.m3 = m3;
    return t;
  }

template <typename T0, typename T1, typename T2, typename T3, typename T4>
  tuple<T0, T1, T2, T3, T4> make_tuple(const T0 & m0, const T1 & m1, const T2 & m2, const T3 & m3, const T4 & m4)
  {
    tuple<T0, T1, T2, T3, T4> t;
    t.m0 = m0;
    t.m1 = m1;
    t.m2 = m2;
    t.m3 = m3;
    t.m4 = m4;
    return t;
  }

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
  tuple<T0, T1, T2, T3, T4, T5> make_tuple(const T0 & m0, const T1 & m1, const T2 & m2, const T3 & m3, const T4 & m4, const T5 & m5)
  {
    tuple<T0, T1, T2, T3, T4, T5> t;
    t.m0 = m0;
    t.m1 = m1;
    t.m2 = m2;
    t.m3 = m3;
    t.m4 = m4;
    t.m5 = m5;
    return t;
  }

} // namespace sps


#endif

