%module elai

%include <mpi4py/mpi4py.i>
%mpi4py_typemap(Comm, MPI_Comm);

%{
#define SWIG_FILE_WITH_INIT
#define ELAI_USE_MPI 1
#define ELAI_USE_MUMPS 1
#define ELAI_USE_PYTHON 1
#define ELAI_FLATNESTED 1
#include <mpi.h>
#include <elai.h>
#include "proxy.hpp"
%}

// Get the NumPy typemaps.
%include "numpy.i"

%init %{
import_array();
%}

%apply (int *DIM1, int **ARGOUTVIEW_ARRAY1) {(int *nelems, int **data)};
%apply (int *DIM1, double **ARGOUTVIEW_ARRAY1) {(int *nelems, double **data)};
%apply (int *IN_ARRAY1, int DIM1) {(int *row, int nrows)};
%apply (int *IN_ARRAY1, int DIM1) {(int *col, int ncols)};
%apply (int DIM1, int *IN_ARRAY1) {(int ncols, int *col)};
%apply (int DIM1, int *IN_ARRAY1) {(int nelems, int *data)};
%apply (int DIM1, double *IN_ARRAY1) {(int nelems, double *data)};
%apply (int DIM1, int DIM2, double *INPLACE_ARRAY2) {(int nrows, int ncols, double *data)};

#define ELAI_USE_MPI 1
#define ELAI_USE_MUMPS 1
#define ELAI_USE_PYTHON 1
#define ELAI_FLATNESTED 1

%include "proxy.hpp"
%include <Elai/sync.hpp>
%include <Elai/ksp.hpp>
%include <Elai/preconditioner.hpp>
%include <Elai/coherence.hpp>
%include <Elai/matrix.hpp>
%include <Elai/vector.hpp>
%include <Elai/bicgstab.hpp>
%include <Elai/bicgsafe.hpp>
%include <Elai/gmres.hpp>
%include <Elai/ilu.hpp>
%include <Elai/subjugator.hpp>
%include <Elai/portal.hpp>
%include <Elai/linear_function.hpp>
%include <Elai/linear_operator.hpp>
%include <Elai/space.hpp>
%include <Elai/family.hpp>
%include <Elai/lu.hpp>
%include <Elai/generator.hpp>
%include <Elai/mumps.hpp>
%include "std_vector.i"

%template(ElaiMatrixDouble) elai::matrix<double>;
%template(ElaiVectorDouble) elai::vector<double>;
%template(ElaiKspDouble) elai::ksp<double>;
%template(ElaiBicgstabDouble) elai::bicgstab<double>;
%template(ElaiBicgsafeDouble) elai::bicgsafe<double>;
%template(ElaiGmresDouble) elai::gmres<double>;
%template(ElaiIluDouble)    elai::ilu<double>;
%template(ElaiGeneratorDouble) elai::generator<double>;
%template(ElaiFunctionDouble) elai::linear_function< elai::Element, elai::Neighbour, double >;
%template(ElaiSubjugatorDouble) elai::subjugator< elai::Element, elai::Neighbour >;
%template(ElaiOperatorDouble) elai::linear_operator< elai::Element, elai::Neighbour, double >;
%template(ElaiSpaceDouble) elai::space<elai::Element>;
%template(ElaiFamilyDouble) elai::family<elai::Element, elai::Neighbour>;
%template(ElaiLuDouble)    elai::lu<double>;
%template(ElaiMumpsDouble) elai::mumps<double>;

%extend elai::coherence {
    %template(ElaiCoherenceDouble) coherence<elai::linear_function< elai::Element, elai::Neighbour, double >>;
};

%extend elai::sync {
    %template(ElaiSyncDouble) sync<elai::linear_function< elai::Element, elai::Neighbour, double >>;
};
%rename(__call__) elai::sync::operator();

%inline %{
    elai::generator<double>::Matrix* ElaiCastToMatrixDouble(elai::matrix<double>* ptr)
    {
        return static_cast<elai::generator<double>::Matrix* >(ptr);
    }

    elai::space<elai::Element> *
        ElaiCastToSpaceGeneratorElementDouble(
            elai::generator<double>::Space* ptr)
    {
        return static_cast<elai::space<elai::Element> *>(ptr);
    }

    elai::family<elai::Element, elai::Neighbour> *
        ElaiCastToFamilyGeneratorElementDouble(
            elai::generator<double>::Family* ptr)
    {
        return static_cast<elai::family<elai::Element, elai::Neighbour> *>(ptr);
    }

    elai::coherence* ElaiCastToCoherence(elai::coherence *ptr)
    {
        return static_cast<elai::coherence *>(ptr);
    }

%}

%extend elai::linear_function< elai::Element, elai::Neighbour, double > {
    inline void setVector(const Element& x, double val)
    {
        (*$self)(x) = val;
    }
}

%extend elai::linear_operator< elai::Element, elai::Neighbour, double > {
    inline void setMatrix(const Element& f, double val)
    {
	(*$self)(f) = val;
    }
}

%extend elai::linear_operator< elai::Element, elai::Neighbour, double > {
    inline void setMatrix(const Element& f, const Element& x, double val)
    {
	(*$self)(f, x) = val;
    }
}

namespace std {
    %template(vectorInt) vector<int>;
};
