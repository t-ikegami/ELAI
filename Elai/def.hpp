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
#ifndef __ELAI_DEF__
#define __ELAI_DEF__

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include "config.hpp"
#include "util.hpp"

#ifdef ELAI_USE_MPI
#include "mpi.h"
#endif

#define ELAI_MAJOR 0.6
#define ELAI_MINOR 2
#define ELAI_VERSION ELAI_MAJOR.ELAIMINOR

#if 201112L <= __STDC_VERSION__
#define ELAI_USE_C11
#endif

namespace elai
{

typedef std::size_t size_t;

template< bool, typename, typename > struct If_;
template< typename T1, typename T2 > struct If_< true, T1, T2 > { typedef T1 type; };
template< typename T1, typename T2 > struct If_< false, T1, T2 > { typedef T2 type; };

template< typename Coef >
Coef conj_( const Coef& v ) { return v; }
template< typename Coef >
std::complex< Coef > conj_( const std::complex< Coef >& v ) { return std::conj( v ); }

#ifdef ELAI_USE_MPI
template< typename T > class mpi_
{
public:
  MPI_Datatype type;
  mpi_( T base ) : type( static_cast< MPI_Datatype >( base ) ) {}
};
template<> class mpi_< bool >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_LOGICAL ) {}
};
template<> class mpi_< int >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_INT ) {}
};
template<> class mpi_< long long >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_LONG ) {}
};
template<> class mpi_< float >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_FLOAT ) {}
};
template<> class mpi_< double >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_DOUBLE ) {}
};
template<> class mpi_< std::complex< float > >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_COMPLEX ) {}
};
template<> class mpi_< std::complex< double > >
{
public:
  MPI_Datatype type;
  mpi_() : type( MPI_DOUBLE_COMPLEX ) {}
};
/*
template<> struct mpi_datatype< int > { static MPI_Datatype type; };
template<> struct mpi_datatype< long long > { static MPI_Datatype type; };
template<> struct mpi_datatype< float > { static MPI_Datatype type; };
template<> struct mpi_datatype< double > { static MPI_Datatype type; };
template<> struct mpi_datatype< std::complex< float > > { static MPI_Datatype type; };
template<> struct mpi_datatype< std::complex< double > > { static MPI_Datatype type; };
MPI_Datatype mpi_datatype< int >::type = MPI_INT;
MPI_Datatype mpi_datatype< long long >::type = MPI_LONG;
MPI_Datatype mpi_datatype< float >::type = MPI_FLOAT;
MPI_Datatype mpi_datatype< double >::type = MPI_DOUBLE;
MPI_Datatype mpi_datatype< std::complex< float > >::type = MPI_COMPLEX;
MPI_Datatype mpi_datatype< std::complex< double > >::type = MPI_DOUBLE_COMPLEX;
*/
#endif//ELAI_USE_MPI

}

#endif//__ELAI_DEF__
