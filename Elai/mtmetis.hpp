/*
 *
 * Elastic Linear Algebra Interface (ELAI)
 *
 * Copyright 2013-2015 H. KOSHIMOTO, AIST
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#ifndef __ELAI_MTMETIS__
#define __ELAI_MTMETIS__

#include <vector>
#include "mtmetis.h"
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"

namespace elai
{

template< class Element, class Neighbour >
class mtmetis
{
public:
  typedef space< Element > Space;
  typedef struct Space::point Point;
  typedef struct Space::const_point CPoint;
  typedef family< Element, Neighbour > Family;

private:
  const Space& base_;
  const Family& tau_;

  int nvtxs_;
  int *xadj_, *adjy_;
  int *perm_; // A'(i) = A(perm_[i])
  int *iperm_;// A(i) = A'(iperm_[i])

public:
  mtmetis( const Space& base, const Family& tau )
    : base_( base ), tau_( tau ), nvtxs_( base_.size() )
  {
    xadj_ = new int[ nvtxs_ + 1 ];

    std::vector< Space > adjs;
    int offset = 0;
    xadj_[ 0 ] = 0;
    for ( typename Space::const_iterator it = base_.begin()
        ; it != base_.end(); ++it
        )
    { // Be care of retreaving diagonals.
      const Space& s = tau_( CPoint( it ).element );
      Space p; p.join( CPoint( it ).element );

      adjs.push_back( s / p );
      xadj_[ CPoint( it ).index + 1 ] = adjs[ offset++ ].size();
    }
    for ( int i = 0; i < nvtxs_; ++i ) xadj_[ i + 1 ] += xadj_[ i ];

    adjy_ = new int[ xadj_[ nvtxs_ ] ];
    offset = 0;
    for ( typename Space::const_iterator it = base_.begin()
        ; it != base_.end(); ++it
        )
    {
      const Space& adj = adjs[ offset++ ];
      int i = CPoint( it ).index;
      int k = 0;

      for ( typename Space::const_iterator jt = adj.begin()
          ; jt != adj.end(); ++jt
          )
      {
        int j = base_.index( CPoint( jt ).element );

        adjy_[ xadj_[ i ] + k ] = j;
        ++k;
      }
    }
    perm_ = new int[ nvtxs_ ];
    iperm_ = new int[ nvtxs_ ];

    reorder();
  }
  ~mtmetis ()
  { delete [] iperm_; delete [] perm_; delete [] adjy_; delete [] xadj_; }

  void reorder() { mtmteis_nd( nvtxs_, xadj_, adjy_, NULL, NULL, perm_ ); }

  Space reordered()
  {
    Space s( base_ );

    for ( typename Space::iterator it = s.begin(); it != s.end(); ++it )
      Point( it ).index = perm_[ Point( it ).index ];

    return s;
  }

  Space inversed( const Space& base ) const
  {
    Space s( base );

    for ( typename Space::iterator it = s.begin(); it != s.end(); ++it )
      Point( it ).index = iperm_[ Point( it ).index ];

    return s;
  }

  const int *perm() const { return perm_; }
  const int *iperm() const { return iperm_; }
};

}

#endif//__ELAI_MTMETIS__
