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
#ifndef __ELAI__
#define __ELAI__

#include "Elai/config.hpp"
#include "Elai/def.hpp"
#include "Elai/coherence.hpp"
#include "Elai/sync.hpp"
#include "Elai/portal.hpp"
#include "Elai/expression.hpp"
#include "Elai/vector.hpp"
#include "Elai/matrix.hpp"
#include "Elai/blas.hpp"
#include "Elai/space.hpp"
#include "Elai/family.hpp"
#include "Elai/clique.hpp"
#include "Elai/subjugator.hpp"
#include "Elai/generator.hpp"
#include "Elai/linear_function.hpp"
#include "Elai/entire_function.hpp"
#include "Elai/linear_operator.hpp"
#include "Elai/entire_operator.hpp"
#include "Elai/preconditioner.hpp"
#include "Elai/fillin.hpp"
#include "Elai/ksp.hpp"
#include "Elai/jacobi.hpp"
#include "Elai/sor.hpp"
#include "Elai/cg.hpp"
#include "Elai/bicgstab.hpp"
#include "Elai/bicgsafe.hpp"
#include "Elai/gmres.hpp"
#include "Elai/jacobi_conditioner.hpp"
#include "Elai/sor_conditioner.hpp"
#include "Elai/ic.hpp"
#include "Elai/ilu.hpp"
#include "Elai/lu.hpp"

#ifdef ELAI_USE_METIS
#include "Elai/metis.hpp"
#endif

#ifdef ELAI_USE_MTMETIS
#include "Elai/mtmetis.hpp"
#endif

#ifdef ELAI_USE_MUMPS
#include "Elai/mumps.hpp"
#endif

#ifdef ELAI_USE_SUPERLU
#include "Elai/superlu.hpp"
#endif

#endif//__ELAI__
