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
#ifndef __ELIA_CONFIG__
#define __ELIA_CONFIG__

/*
 * User Configurations
 */
//#define ELAI_DEBUG
//#define ELAI_PROFILE
//#define ELAI_USE_OPENMP
//#define ELAI_USE_MPI
//#define ELAI_USE_MPI3
//#define ELAI_USE_MUMPS
//#define ELAI_USE_SUPERLU
//#define ELAI_USE_METIS
//#define ELAI_USE_MTMETIS

#define ELAI_USE_METIS

#if defined( ELAI_USE_MPI3 ) & !defined( ELAI_USE_MPI )
#define ELAI_USE_MPI
#endif

#ifdef ELAI_DEBUG
#define NDEBUG
#endif//ELAI_DEBUG

namespace elai
{

// Alignment-size depends on your architecture.
// Default: Intel64 on Haswell
const size_t DATA_ALIGNMENT = 16;
const size_t CACHE_LINE = 64;

}

#endif//__ELIA_CONFIG__
