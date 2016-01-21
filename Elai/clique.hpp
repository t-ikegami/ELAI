/*
 *
 * Enexss Linear Algebra Interface (ELAI)
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
#ifndef __ELAI_CLIQUE__
#define __ELAI_CLIQUE__

#include <vector>
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace elai
{

template< class Element, class Neighbour >
class clique
{
  const space< Element >& s_;
  const family< Element, Neighbour >& tau_;

  int n_, nnz_;
  int *xadj_;
  int *adjy_;

  void terminate()
  {
    if ( xadj_ != NULL ) { delete [] xadj_; xadj_ = NULL; }
    if ( adjy_ != NULL ) { delete [] adjy_; adjy_ = NULL; }
  }

public:
  clique( const space< Element >& s, const family< Element, Neighbour >& tau )
    : s_( s ), tau_( tau ), n_( s.size() ), nnz_( 0 ), xadj_( NULL ), adjy_( NULL )
  {}
  ~clique()
  { terminate(); }

  void fill( int lv )
  {
    std::vector< space< Element > > adjs;
    int offset;

    terminate();
    xadj_ = new int[ n_ + 1 ];
    xadj_[ 0 ] = 0;
    offset = 0;
    for ( typename space< Element >::const_iterator it = s_.begin()
        ; it != s_.end(); ++it
        )
    {
      adjs.push_back( tau_( Point( it ).element, lv ) );
      xadj_[ Point( it ).index + 1 ] = adjs[ offset++ ].size();
    }
    for ( int i = 0; i < n_; ++i ) xadj_[ i + 1 ] += xadj_[ i ];
    nnz_ = xadj_[ n_ ];

    adjy_ = new int[ nnz_ ];
    offset = 0;
    for ( typename space< Element >::const_iterator it = s_.begin()
        ; it != s_.end(); ++it
        )
    {
      space< Element >& adj = adjs[ offset++ ];
      int i = Point( it ).index;
      int k = 0;

      for ( typename space< Element >::const_iterator jt = adj.begin()
          ; jt != adj.end(); ++jt
          )
      {
        int j = s_.index( Point( jt ).element );

        adjy_[ xadj_[ i ] + k ] = j;
        ++k;
      }
    }
  }

  int n() const { return n_; }
  int nnz() const { return nnz_; }
  int *xadj() const { return xadj_; }
  int *adjy() const { return adjy_; }
};

}

#endif//__ELAI_CLIQUE__
