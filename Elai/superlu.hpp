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
#ifndef __ELAI_SUPERLU__
#define __ELAI_SUPERLU__

#include "def.hpp"
#include "lu.hpp"

extern "C"
{
#include "superlu_dist.h"
}

namespace elai
{

template < class T > struct mumps_impl;
template <> struct mumps_impl < float >
{
  typedef SMUMPS_STRUC_C type;
  static void call( type *mumps ) { smumps_c( mumps ); }
};
template <> struct mumps_impl < double >
{
  typedef DMUMPS_STRUC_C type;
  static void call( type *mumps ) { dmumps_c( mumps ); }
};
template <> struct mumps_impl < std::complex < float > >
{
  typedef CMUMPS_STRUC_C type;
  static void call( type *mumps ) { cmumps_c( mumps ); }
};
template <> struct mumps_impl < std::complex < double > >
{
  typedef ZMUMPS_STRUC_C type;
  static void call( type *mumps ) { zmumps_c( mumps ); }
};

template < class Coef >
class mumps : public lu< Coef >
{
  ELAI_USE_LU;

  typedef mumps_impl < Coef > Mumps;
  typedef typename Mumps::type MumpsType;

  MumpsType mumps_;

  bool factor_()
  {
    mumps_.job = 2;

    if ( rank_ == 0 ) mumps_.a = A_.val();
    //else mumps_.a_loc = A_.val();
    Mumps::call( &mumps_ );

    mem_ = mumps_.INFO( 16 ); // Elapsed USAGE

    return 0 <= mumps_.INFO( 1 );
  }

  bool solve_( const vector< Coef >& b, vector< Coef >& x )
  {
    mumps_.job = 3;

    if ( rank_ == 0 )
    {
      for ( int i = 0 ; i < x.m(); ++i ) x( i ) = b( i );
      mumps_.rhs = x.val();
    }
    Mumps::call( &mumps_ );

    return 0 <= mumps_.INFO( 1 );
  }

public:
  mumps( matrix< Coef >& A, MPI_Comm comm )
    : lu< Coef >( A, comm ), mumps_()
  {
    mumps_.job = JOB_INIT;
    mumps_.par = 1;
    mumps_.sym = 0;
    mumps_.comm_fortran = MPI_Comm_c2f( comm_ );
    Mumps::call( &mumps_ );

    // Analyse
    mumps_.job = 1;
    mumps_.CNTL( 1 ) = 0e0;

    // Log level
    //mumps_.ICNTL( 1 ) = 6;
    //mumps_.ICNTL( 2 ) = 0;
    //mumps_.ICNTL( 3 ) = 6;
#ifdef ELAI_DEBUG
    mumps_.ICNTL( 4 ) = 4;
#else
    mumps_.ICNTL( 4 ) = 1;
#endif

    mumps_.ICNTL( 10 ) = 2; // #of Iterative Refinement
    mumps_.ICNTL( 11 ) = 0; // Error Analysis
    //if ( size_ == 1 )
    if ( 1 )
    { // SINGLE
      mumps_.ICNTL( 5 ) = 0; // Assembled Format
      mumps_.ICNTL( 6 ) = 7; // Unsymmetric Permutation
      //mumps_.ICNTL( 7 ) = 1; // USER
      //mumps_.ICNTL( 7 ) = 2; // AMF
      mumps_.ICNTL( 7 ) = 5; // METIS
      mumps_.ICNTL( 18 ) = 0; // Storage Centralized
      mumps_.ICNTL( 21 ) = 0; // Solution Centralized
      if ( rank_ == 0 )
      {
        mumps_.n = A_.m();
        mumps_.nz = A_.nnz();
        mumps_.irn = new int[ mumps_.nz ];
        for ( int i = 0, off = 0; i < A_.m(); ++i )
          for ( int k = A_.ind( i ); k < A_.ind( i + 1 ); ++k )
            mumps_.irn[ off++ ] = i + 1;
        mumps_.jcn = new int[ mumps_.nz ];
        for ( int i = 0, off = 0; i < A_.m(); ++i )
          for ( int k = A_.ind( i ); k < A_.ind( i + 1 ); ++k )
            mumps_.jcn[ off++ ] = A_.col( k ) + 1;
        //mumps_.perm_in = perm; // Ordering
      }
    }
    else
    { // MULTI for next-release
      mumps_.ICNTL( 5 ) = 0; // Assembled Format
      mumps_.ICNTL( 6 ) = 0; // Unsymmetric Permutation
      //mumps_.ICNTL( 7 ) = 1; // USER
      //mumps_.ICNTL( 7 ) = 2; // AMF
      mumps_.ICNTL( 7 ) = 5; // METIS
      mumps_.ICNTL( 18 ) = 3; // Storage Distributed
      mumps_.ICNTL( 20 ) = 0;
      mumps_.ICNTL( 21 ) = 0; // Solution Centralized
      if ( rank_ == 0 )
      {
        mumps_.n = A_.m();
        mumps_.nz = A_.nnz();
        mumps_.irn = new int[ mumps_.nz ];
        for ( int i = 0, off = 0; i < A_.m(); ++i )
          for ( int k = A_.ind( i ); k < A_.ind( i + 1 ); ++k )
            mumps_.irn[ off++ ] = i + 1;
        mumps_.jcn = new int[ mumps_.nz ];
        for ( int i = 0, off = 0; i < A_.m(); ++i )
          for ( int k = A_.ind( i ); k < A_.ind( i + 1 ); ++k )
            mumps_.jcn[ off++ ] = A_.col( k ) + 1;
        //mumps_.perm_in = perm; // Ordering
        //mumps_.nz_loc = ?.nnz();
        //mumps_.irn_loc = ?.ind();
        //mumps_.jcn_loc = ?.col();
        //mumps_.mapping = mapping; // Affiliation Ranks
      }
      else
      {
        // TODO
        mumps_.nz_loc = A_.nnz();
        mumps_.irn = NULL;
        mumps_.irn_loc = A_.ind();
        mumps_.jcn = NULL;
        mumps_.jcn_loc = A_.col();
      }
    }
    mumps_.ICNTL( 14 ) = 30; // Work space ( default: 20 )
    mumps_.ICNTL( 24 ) = 1; // Pivot On
    mumps_.ICNTL( 28 ) = 1; // Sequencial ordering
    mumps_.CNTL( 1 ) = 1e-05; // Numerical Pivot Quality ( default: 0.01 )
    mumps_.CNTL( 3 ) = 1e-17; // Dynamic Pivot Threshold
    mumps_.CNTL( 4 ) = 1e-20; // Static Pivot Threshold
    mumps_.CNTL( 5 ) = 1e+18; // Null diagonal scale
    Mumps::call( &mumps_ );

    if ( 1 < size_ )
    {
      int mems[ size_ ];

      MPI_Gather( &mumps_.INFO( 15 ), 1, MPI_INT, mems, 1, MPI_INT, 0, comm_ );
      for ( int i = 0; i < size_; ++i ) mem_ += mems[ i ]; // Approximated USAGE
    }
    else mem_ = mumps_.INFO( 15 );
  }
  ~mumps()
  {
    if ( mumps_.jcn != NULL ) { delete [] mumps_.jcn; mumps_.jcn = NULL; }
    if ( mumps_.irn != NULL ) { delete [] mumps_.irn; mumps_.irn = NULL; }
  }
};

}

#undef JOB_INIT
#undef JOB_END
#undef INFO(I)
#undef RINFO(I)
#undef CNTL(I)
#undef ICNTL(I)

#endif//__ELAI_SUPERLU__
