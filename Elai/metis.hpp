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
#ifndef __ELAI_METIS__
#define __ELAI_METIS__

#include <vector>
#include "metis.h"
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"

namespace elai
{

template< class Element, class Neighbour >
class metis
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
  int *options_;
  /*
   * METIS_OPTION_PTYPE - METIS_PTYPE_RB, METIS_PTYPE_KWAY
   * METIS_OPTION_OBJTYP - METIS_OBJTYPE_CUT, METIS_OBJTYPE_VOL
   * METIS_OPTION_CTYPE - METIS_CTYPE_RM, METIS_CTYPE_SHEM
   * METIS_OPTION_IPTYPE - METIS_IPTYPE_GROW, METIS_IPTYPE_RANDOM, METIS_IPTYPE_EDGE, METIS_IPTYPE_NODE
   * METIS_OPTION_RTYPE - METIS_RTYPE_FM, METIS_RTYPE_GREEDY, METIS_RTYPE_SEP2SIDED, METIS_RTYPE_SEP1SIDED
   * METIS_OPTION_NCUTS - The number of partitioning candidates (Default: 1)
   * METIS_OPTION_NSEPS - The number of separator candidates (Default: 1)
   * METIS_OPTION_NUMBERING - The origin of index (C: 0, Fortran: 1)
   * METIS_OPTION_NITER - The number of refinement iterations (Default: 10)
   * METIS_OPTION_SEED - The random generator
   * METIS_OPTION_MINCONN - Minimize the maximum degree (DONT: 0, DO: 1)
   * METIS_OPTION_NO2HOP - 2-hop matching (DO: 0, DONT: 1)
   * METIS_OPTION_CONTIG - Try to produce contiguous partitions (DONT: 0, DO: 1)
   * METIS_OPTION_COMPRESS - Try to compress the graph (DONT: 0, DO: 1)
   * METIS_OPTION_CCORDER - Identify connected components (DONT: 0, DO: 1)
   * METIS_OPTION_PFACTOR - Remove vertices whose degree is less than 0.1*x*ave (Default: 0)
   * METIS_OPTION_UFACTOR - The allowable ratio of load imbalance (1+x)/1000 (Default: 30)
   * METID_OPTION_DGBLVL
   */

public:
  metis( const Space& base, const Family& tau )
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

    options_ = new int[ METIS_NOPTIONS ];
    METIS_SetDefaultOptions( options_ );
    //options_[ METIS_OPTION_PTYPE ] = METIS_PTYPE_RB;
    //options_[ METIS_OPTION_OBJTYPE ] = METIS_OBJTYPE_CUT;
    //options_[ METIS_OPTION_CTYPE ] = METIS_CTYPE_SHEM;
    //options_[ METIS_OPTION_IPTYPE ] = METIS_IPTYPE_EDGE;
    //options_[ METIS_OPTION_RTYPE ] = METIS_RTYPE_SEP2SIDED;
    //options_[ METIS_OPTION_NUMBERING ] = 0;
    //options_[ METIS_OPTION_NCUTS ] = 1;
    //options_[ METIS_OPTION_NSEPS ] = 1;
    //options_[ METIS_OPTION_NITER ] = 3;

    reorder();
  }
  ~metis ()
  { delete [] options_; delete [] iperm_; delete [] perm_; delete [] adjy_; delete [] xadj_; }

  int& option( int opt ) { return options_[ opt ]; }
  const int& option( int opt ) const { return options_[ opt ]; }
  void reorder() { METIS_NodeND( &nvtxs_, xadj_, adjy_, NULL, options_, perm_, iperm_ ); }

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

#endif//__ELAI_METIS__
