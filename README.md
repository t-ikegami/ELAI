# Elastic Linear Algebra Interface (ELAI)

## Introduction

ELAI is a C++ template library for the parallel sparse linear solver.
It was originally designed for the linear problems in semiconductor
device simulation.  With ELAI library, you can easily take physical
properties of the target problem into account in solving the linear
problems: ordering and partitioning can be performed considering
physical constraints.

ELAI includes but not limited to the following functionalities.
- merging matrices decomposed via domains.
- LU factoring with distributed matrix.
- additive schwarz preconditioner for decomposed matrices.
- structural ordering on coupled systems.
- unsymmetric permutations for maximize the product on diagonals.
- column-wise scalings.

Following libraries are required to use ELAI:
- MPI
  (We have tested with OpenMPI 1.6.5)
- [MUMPS 5.0.1](http://mumps.enseeiht.fr)
  BLAS, LAPACK, and ScaLAPACK are also required.
- [METIS 5.1.0](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)

## License
Copyright 2013-2016 H. KOSHIMOTO, AIST

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
