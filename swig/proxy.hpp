#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Elai/vector.hpp"
#include "Elai/matrix.hpp"
#include "Elai/bicgstab.hpp"

volatile bool debug = true;

void
DebugStop(void)
{
    while (debug) {
        sleep(1);
    }
}

std::vector< int >*
NewStdVector(int nelems, int *data)
{
    assert(nelems > 0);
    assert(data != NULL);
    std::vector <int>* vec = new std::vector<int>(nelems);
    for (int i = 0; i < nelems; i++) {
        (*vec)[i] = data[i];
    }
    return vec;
}

#if 0
// This function is not necessary for Python. Because, Python will call
// destructor when instance is not necessary.
void
DeleteStdVector(std::vector< int >* vec)
{
    assert(vec != NULL);
    delete vec;
}
#endif // 0

std::istream *ProxyIfstream(const char *fname)
{
    return new std::ifstream(fname);
}

void
ProxyResult(
    elai::matrix< double >& A,
    elai::vector< double >& x,
    elai::vector< double >& b)
{
    elai::vector< double> v = A * x - b;
    double r = v * v, r0;

    r = sqrt(r);
    std::cout << "||A x - b|| = " << std::setprecision( 15 ) << r << std::endl;
    r0 = x * x; r0 = sqrt( r0 );
    r /= r0;
    std::cout << "||A x - b||/||x|| = " << std::setprecision( 15 ) << r << std::endl;
}

void
mul(
    elai::matrix< double >& mat,
    elai::vector< double >& vec,
    elai::vector< double >& result)
{
    result = mat * vec;
}

double
mul(
    elai::vector<double >& src,
    elai::vector<double >& dst)
{
    return src * dst;
}

void
minus(
    elai::vector<double >& src,
    elai::vector<double >& dst,
    elai::vector<double >& result)
{
    result = src - dst;
}
