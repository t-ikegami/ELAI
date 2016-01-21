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
#ifndef __ELAI_ENTIRE_FUNCTION__
#define __ELAI_ENTIRE_FUNCTION__

#include <vector>
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"
#include "linear_function.hpp"

namespace elai
{

template< class Element, class Neighbour, class Range >
class entire_function
{
  typedef space< Element > Space;
  typedef family< Element, Neighbour > Family;
  typedef linear_function< Element, Neighbour, Range > Function;

  std::vector< const Function * > handles_;

public:
  entire_function() {}
  ~entire_function() {}

  entire_function< Element, Neighbour, Range >& join( const Function& func )
  {
    handles_.push_back( &func );

    return *this;
  }

  void purge()
  {
    handles_.clear();
  }

  linear_function< Element, Neighbour, Range >
  operator()()
  {
    Space domain;
    Range *v;

    for ( int i = 0; i < handles_.size(); ++i )
      domain |= handles_[ i ]->dom();

    v = new Range[ domain.size() ];
    for ( int i = 0; i < handles_.size(); ++i )
    {
      const Function& f = *handles_[ i ];
      const Space& dom = f.dom();

      for ( typename Space::const_iterator it = dom.begin()
          ; it != dom.end(); ++it
          )
      {
        typename Space::const_point pt( it );

        if ( !dom.govern( pt.element ) ) continue;
        int i = domain.index( pt.element );

        v[ i ] = f( pt.element );
      }
    }
    linear_function< Element, Neighbour, Range > func( domain, v );
    delete [] v;

    return func;
  }
};

}

#endif//__ELAI_ENTIRE_FUNCTION__
