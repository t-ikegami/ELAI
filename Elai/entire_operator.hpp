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
#ifndef __ELAI_ENTIRE_OPERATOR__
#define __ELAI_ENTIRE_OPERATOR__

#include <vector>
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"
#include "matrix.hpp"
#include "linear_operator.hpp"

namespace elai
{

template< class Element, class Neighbour, class Range >
class entire_operator
{
  typedef space< Element > Space;
  typedef family< Element, Neighbour > Family;
  typedef linear_operator< Element, Neighbour, Range > Operator;
  typedef matrix< Range > Matrix;

  std::vector< const Operator * > operators_;

public:
  entire_operator() {}
  ~entire_operator() {}

  entire_operator< Element, Neighbour, Range >& join( const Operator& op )
  {
    operators_.push_back( &op );

    return *this;
  }

  void purge()
  {
    operators_.clear();
  }

  linear_operator< Element, Neighbour, Range >
  operator()( bool is_symmetric = false )
  {
    Space domain;
    Space range;
    Family graph;

    for ( int i = 0; i < operators_.size(); ++i )
    {
      const Operator& op = *operators_[ i ];
      Family topo = op.topo();

      for ( typename Space::const_marginal_iterator it = op.ran().internal_begin()
          ; it != op.ran().internal_end(); ++it
          )
      {
        topo.flip
          ( typename Space::const_internal_point( it ).internal
          , typename Space::const_internal_point( it ).external
          , is_symmetric
          );
      }
      domain |= op.dom();
      range |= op.ran();
      graph |= topo;
    }
    // THE LAST REORDERING CHANCE

    std::vector< Space > adjs;
    int *ind = new int[ range.size() + 1 ], m = range.size(), off;

    off = 0;
    ind[ 0 ] = 0;
    for ( typename Space::const_iterator it = range.begin(); it != range.end(); ++it )
    {
      const Element& e = typename Space::const_point( it ).element;
      int index = typename Space::const_point( it ).index;

      adjs.push_back( graph( e ) );
      ind[ index + 1 ] = adjs[ off++ ].size();
    }
    for ( int i = 0; i < m; ++i ) ind[ i + 1 ] += ind[ i ];

    int nnz = ind[ m ];
    int *col = new int[ nnz ];
    Range *c = new Range[ nnz ];

    off = 0;
    for ( typename Space::const_iterator it = range.begin(); it != range.end(); ++it )
    {
      const Space& adj = adjs[ off++ ];
      const Element& e = typename Space::const_point( it ).element;
      int i = typename Space::const_point( it ).index;
      int k = 0;

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& d = typename Space::const_point( jt ).element;
        int j = domain.index( d );

        col[ ind[ i ] + k ] = j;
        ++k;
      }
    }

    for ( int n = 0; n < operators_.size(); ++n )
    {
      const Operator& op = *operators_[ n ];
      const Space& ran = op.ran();
      const Space& dom = op.dom();
      const Family& topo = op.topo();
      const Matrix& A = op.action();

      for ( typename Space::const_iterator it = ran.begin()
          ; it != ran.end(); ++it
          )
      {
        const Element& x = typename Space::const_point( it ).element;
        const Space& adj = topo( x );
        int I = typename Space::const_point( it ).index;

        if ( ran.govern( x ) )
        {
          int i = range.index( x );

          for ( typename Space::const_iterator jt = adj.begin()
              ; jt != adj.end(); ++jt
              )
          {
            const Element& y = typename Space::const_point( jt ).element;
            int J = dom.index( y );
            int j;

            if ( dom.govern( y ) ) j = domain.index( y );
            else j = domain.index( dom.external_element( y ) );

            for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
              if ( col[ k ] == j )
              {
                c[ k ] = A( I, J );
                break;
              }
          }
        }
      }
    }

    Matrix A( range.size(), domain.size(), nnz, ind, col, c );

    delete [] c;
    delete [] col;
    delete [] ind;

    return linear_operator< Element, Neighbour, Range >( range, domain, graph, A );
  }
};

}

#endif//__ELAI_ENTIRE_OPERATOR__
