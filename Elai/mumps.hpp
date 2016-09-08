/*
 *
 * Elastic Linear Algebra Interface (ELAI)
 *
 * Copyright 2013-2016 H. KOSHIMOTO, AIST
 * Copyright 2016 T. IKEGAMI, AIST
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
#ifndef __ELAI_MUMPS__
#define __ELAI_MUMPS__

#include "def.hpp"
#include "lu.hpp"

#define JOB_INIT	-1
#define JOB_END		-2
#define JOB_ANALYSE	 1
#define JOB_FACTORIZE	 2
#define JOB_SOLVE	 3
#define INFO(I)   info[(I)-1]		/* C <-> FORTRAN index conversion */
#define INFOG(I)  infog[(I)-1]
#define RINFO(I)  rinfo[(I)-1]
#define RINFOG(I) rinfog[(I)-1]
#define CNTL(I)   cntl[(I)-1]
#define ICNTL(I)  icntl[(I)-1]

extern "C"
{
#include "smumps_c.h"
#include "dmumps_c.h"
#include "cmumps_c.h"
#include "zmumps_c.h"
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

struct mumps_options {
public:
  int    par, sym, debug, iter, stat, work, bcast;
  double pivot_quality, pivot_dynamic, pivot_static, pivot_fixation;

  mumps_options() :
    par(1),		// master rank is involved in the calculation
    sym(0),		// 0:asymmetric, 1:symmetric positive definite, 2:symmetric
    debug(1),		// 0:no output, 1:error, 2:warnings&statistics, 3:diagnostics, 4:in/out params
    iter(2),		// number of iterative refinement
    stat(0),		// perform error analysis: 0:none, 1:expensive, 2:moderate
    work(30),		// percentage increase of estimated working space (default: 20)
    bcast(0),		// whether or not the solution vector is bcasted over all ranks
    pivot_quality(1e-05),	// Numerical Pivot Quality (default: 0.01)
    pivot_dynamic(1e-17),	// Dynamic Pivot Threshold (default: 0.0)
    pivot_static(1e-20),	// Static Pivot Threshold  (default: -1.0)
    pivot_fixation(1e+18)	// Null diagonal scale (default: 0.0)
  {}
};

template < class Coef >
class mumps : public lu< Coef >
{
  ELAI_USE_LU;

  typedef mumps_impl < Coef > Mumps;
  typedef typename Mumps::type MumpsType;

  MumpsType mumps_;
  int bcast_;

  bool factor_()
  {
    mumps_.job = JOB_FACTORIZE;

    if ( rank_ == 0 ) mumps_.a = A_.val();
    //else mumps_.a_loc = A_.val();
    Mumps::call( &mumps_ );

    mem_ = mumps_.INFOG( 19 );		// sum of INFO(16); total memory usage in MB

    return 0 <= mumps_.INFOG( 1 );
  }

  bool solve_( const vector< Coef >& b, vector< Coef >& x )
  {
    mumps_.job = JOB_SOLVE;

    if ( rank_ == 0 )
    {
      for ( int i = 0 ; i < x.m(); ++i ) x( i ) = b( i );
      mumps_.rhs = x.val();
    }
    Mumps::call( &mumps_ );

    if ( 0 <= mumps_.INFOG(1) && bcast_ ) {
      MPI_Comm comm = MPI_Comm_f2c(mumps_.comm_fortran);
      MPI_Bcast( x.val(), x.m(), mpi_<Coef>().type, 0, comm );
    }
    
    return 0 <= mumps_.INFOG( 1 );
  }

public:
  mumps( matrix< Coef >& A, MPI_Comm comm, mumps_options opt = mumps_options() )
    : lu< Coef >( A, comm ), mumps_()
  {
    // Initialize
    bcast_ = opt.bcast;
    mumps_.job = JOB_INIT;
    mumps_.par = opt.par;
    mumps_.sym = opt.sym;
    mumps_.comm_fortran = MPI_Comm_c2f( comm_ );
    Mumps::call( &mumps_ );

    // copy matrix
    if ( rank_ == 0 ) {
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
      //mumps_.perm_in = perm;	// user ordering in case ICNTL(7) = 1
    }

    // Analyse
    mumps_.job = JOB_ANALYSE;
    mumps_.ICNTL(  4 ) = opt.debug;	// log level
    mumps_.ICNTL(  5 ) = 0;		// Assembled Format
    mumps_.ICNTL(  6 ) = 7;		// permutation to avoid diagonal zero (automatic choice)
    mumps_.ICNTL(  7 ) = 5;		// ordering (METIS, in case ICNTL(28) = 1)
    mumps_.ICNTL( 10 ) = opt.iter;	// #of Iterative Refinement
    mumps_.ICNTL( 11 ) = opt.stat;	// Error Analysis
    mumps_.ICNTL( 14 ) = opt.work;	// Work space ( default: 20 )
    mumps_.ICNTL( 24 ) = 1;		// Pivot On
    mumps_.ICNTL( 28 ) = 1;		// Sequencial ordering
    mumps_.ICNTL( 29 ) = 2;		// ParMetis ordering (in case ICNTL(28) = 2) 

    mumps_.CNTL( 1 ) = opt.pivot_quality;	// Numerical Pivot Quality
    mumps_.CNTL( 3 ) = opt.pivot_dynamic;	// Dynamic Pivot Threshold
    mumps_.CNTL( 4 ) = opt.pivot_static;	// Static Pivot Threshold
    mumps_.CNTL( 5 ) = opt.pivot_fixation;	// Null diagonal scale

    //if ( size_ == 1 )
    if ( 1 )
    { // SINGLE
      mumps_.ICNTL( 18 ) = 0;		// Storage Centralized
      mumps_.ICNTL( 20 ) = 0;		// dense RHS
      mumps_.ICNTL( 21 ) = 0;		// Solution Centralized
    }
#if 0
    else	// TODO
    { // MULTI (for future-release)
      mumps_.ICNTL( 18 ) = 3;		// Storage Distributed (user defined directly)
      mumps_.ICNTL( 20 ) = 0;		// dense RHS
      mumps_.ICNTL( 21 ) = 0;		// Solution Centralized
      if ( rank_ == 0 )
      {
        mumps_.nz_loc = ?.nnz();
        mumps_.irn_loc = ?.ind();
        mumps_.jcn_loc = ?.col();
        mumps_.mapping = mapping; // Affiliation Ranks
      }
      else
      {
        mumps_.nz_loc = A_.nnz();
        mumps_.irn = NULL;
        mumps_.irn_loc = A_.ind();
        mumps_.jcn = NULL;
        mumps_.jcn_loc = A_.col();
      }
    }
#endif
    Mumps::call( &mumps_ );

    mem_ = mumps_.INFOG(17);		// sum of INFO(15); estimated memory usage in MB
  }
  
  ~mumps()
  {
    if ( mumps_.jcn != NULL ) { delete [] mumps_.jcn; mumps_.jcn = NULL; }
    if ( mumps_.irn != NULL ) { delete [] mumps_.irn; mumps_.irn = NULL; }
    mumps_.job = JOB_END;
    Mumps::call( &mumps_ );
  }

  inline int info (int idx) { return mumps_.INFO(idx); }
  inline int infog(int idx) { return mumps_.INFOG(idx); }
  inline double rinfo (int idx) { return mumps_.RINFO(idx); }
  inline double rinfog(int idx) { return mumps_.RINFOG(idx); }
  
};

}

#undef JOB_INIT
#undef JOB_END
#undef JOB_ANALYSE
#undef JOB_FACTORIZE
#undef JOB_SOLVE
#undef INFO
#undef INFOG
#undef RINFO
#undef RINFOG
#undef CNTL
#undef ICNTL

#endif//__ELAI_MUMPS__
