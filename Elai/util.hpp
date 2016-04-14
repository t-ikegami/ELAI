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
#ifndef __ELAI_UTIL__
#define __ELAI_UTIL__

#include "config.hpp"

extern "C"
{
#include <sys/time.h>
}

#ifdef ELAI_DEBUG

#define ELAI_CHECK( stmt )

#else

#define ELAI_CHECK( stmt )

#endif//ELAI_DEBUG

#ifdef ELAI_PROFILE

#define ELAI_PROF_BEG( acc ) \
  do {                       \
    time_monitor< double > t

#define ELAI_PROF_END( acc ) \
    ( acc ) += t();          \
  } while ( 0 )

#else

#define ELAI_PROF_BEG( acc )
#define ELAI_PROF_END( acc )

#endif//ELAI_PROFILE

namespace elai
{

template< typename Coef >
class time_monitor
{
  Coef origin_;

public:
  time_monitor() : origin_()
  {
    struct timeval tv;

    gettimeofday( &tv, NULL );
    origin_ = static_cast< Coef >( tv.tv_usec / 1000000.0 + tv.tv_sec );
  }

  Coef operator()() const
  {
    struct timeval tv;
    Coef elapsed;

    gettimeofday( &tv, NULL );
    elapsed = static_cast< Coef >( tv.tv_usec / 1000000.0 + tv.tv_sec - origin_ );

    return elapsed;
  }
};

}

#endif//__ELAI_UTIL__
