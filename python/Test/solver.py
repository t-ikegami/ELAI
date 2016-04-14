#!/usr/bin/env python
import argparse
import numpy as np
from mpi4py import MPI
import elai

# Parse arguments.
parser = argparse.ArgumentParser(description = 'Test ELAI solver')
parser.add_argument('afile', nargs = '?', help = 'matrix file')
parser.add_argument('bfile', nargs = '?', help = 'vector file')
parser.add_argument('--method', choices = [ 'mumps', 'gmres', 'bicgsafe', 'bicgstab' ], default = 'mumps',
                    help = 'solver method')
parser.add_argument('--scaled',       action = 'store_true', help = 'enable both-side scaling')
parser.add_argument('--precondition', action = 'store_true', help = 'enalbe preconditioner')
parser.add_argument('--cthres', type = float, default = 1e-08,
                    help = 'relative residual toerance in KSP solvers')
parser.add_argument('--sthres', type = float, default = 1e-05,
                    help = 'threshold for both-side-scaling (diagonal becomes 1 + sthres after the scaling)')
parser.add_argument('--fthres', type = float, default = 1e-08,
                    help = 'threshold for ilu preconditioner')
parser.add_argument('--flevel', type = int, default = 1,
                    help = 'fill-in level for ilu preconditioner')

arg = parser.parse_args()

def solve(A, x, b, coherent = None):
    rnorm = 1.0
    cnorm = 1.0
    ratio = 1.0
    prec = None
    if arg.scaled:
        A.normalize(arg.sthres, coherent)
        b.scale(A.scaleRow())
        x.unscale(A.scaleCol())
        rnorm = A.scaleRowNorm()
        cnorm = A.scaleColNorm()
        ratio = A.scaleRatio()

    if arg.precondition:
        # prec = elai.ElaiIluDouble(A, arg.flevel, arg.fthres, False)
        prec = elai.ElaiIluDouble(A, arg.flevel, arg.fthres)
        prec.factor()

    svr = elai.ElaiBicgstabDouble(A, b, prec, coherent)		if arg.method == 'bicgstab' else \
          elai.ElaiBicgsafeDouble(A, b, prec, coherent)		if arg.method == 'bicgsafe' else \
          elai.ElaiGmresDouble   (A, b, prec, coherent)		if arg.method == 'gmres'    else \
          None
    assert(svr != None)

    svr.iter_max(imax)
    svr.rel_thres(arg.cthres)
    svr.solve(x)

    if arg.scaled:
        x.scale(A.scaleCol())
        b.unscale(A.scaleRow())
        A.unnormalize()

    if MPI.COMM_WORLD.Get_rank() == 0:
        print('SCAL: ' + str(ratio) + ' (' + str(rnorm) + ', ' + str(cnorm) + ')')

def printResult(A, x, b):
    res = elai.ElaiVectorDouble(x)
    elai.mul(A, x, res)
    elai.minus(res, b, res)
    r = elai.mul(res, res)
    r = np.sqrt(r)
    print('||A x - b|| = ' + str(r))
    r0 = elai.mul(x, x)
    r0 = np.sqrt(r0)
    r /= r0
    print('||A x - b||/||x|| = ' + str(r))

def solve_dist(A, U, V, ranks):
    myrank = MPI.COMM_WORLD.Get_rank()
    loc = elai.ElaiSubjugatorDouble(A.dom(), A.topo(), ranks)
    base = elai.ElaiSpaceDouble(loc(myrank))
    topo = elai.ElaiFamilyDouble(loc(myrank, base))
    a = elai.ElaiOperatorDouble(base, topo)
    u = elai.ElaiFunctionDouble(base)
    v = elai.ElaiFunctionDouble(base)
    coherent = elai.ElaiCoherenceDouble(u, MPI.COMM_WORLD)
    sync = elai.ElaiSyncDouble(U, MPI.COMM_WORLD)

    a.reflectIn(A, loc)
    u.reflectIn(U, loc)
    v.reflectIn(V, loc)
    solve(a.action(), u.ran(), v.ran(), coherent)
    U.clear(0.0)
    U.reflect(u, loc)
    sync()

def run_dist(A, x, b):
    gen = elai.ElaiGeneratorDouble(A)
    base = gen.space()
    topo = gen.family()
    a = elai.ElaiOperatorDouble(base, base, topo, A)
    u = elai.ElaiFunctionDouble(base, x)
    v = elai.ElaiFunctionDouble(base, b)
    ranks = elai.NewStdVector(range(MPI.COMM_WORLD.Get_size()))
    solve_dist(a, u, v, ranks)
    if MPI.COMM_WORLD.Get_rank() == 0:
        A = a.action()
        x = u.ran()
        b = v.ran()
        printResult(A, x, b)

def run(A, x, b):
    solve(A, x, b)
    printResult(A, x, b)

def direct(A, x, b):
    lu = elai.ElaiMumpsDouble(A, MPI.COMM_WORLD)
    lu.factor()
    lu.solve(b, x)
    if MPI.COMM_WORLD.Get_rank() == 0:
        print('USE ' + str(MPI.COMM_WORLD.Get_size()) + '-PROC.')
        printResult(A, x, b)

def demo():
    nnd = 16
    color = 0
    mesh = elai.ElaiSpaceDouble()
    topo = elai.ElaiFamilyDouble()
    for i in range(nnd):
        n = elai.Element(i, color)
        mesh.join(n)
        lnk = elai.Neighbour(n)
        if 0 <= (i - 4):
            lnk.join(elai.Element(i - 4, color))
        if 0 <= (i - 1):
            lnk.join(elai.Element(i - 1, color))
        if (i + 1 < nnd):
            lnk.join(elai.Element(i + 1, color))
        if (i + 4 < nnd):
            lnk.join(elai.Element(i + 4, color))
        topo.join(lnk)

    problem = elai.ElaiOperatorDouble(mesh, topo)
    solution = elai.ElaiFunctionDouble(mesh)
    inhomogeneous = elai.ElaiFunctionDouble(mesh)

    for i in range(nnd):
        problem.setMatrix(elai.Element(i, color), -4)
        if 0 <= (i - 4):
            n1 = elai.Element(i - 4, color)
            n2 = elai.Element(i, color)
            problem.setMatrix(n1, n2, 1)
        if 0 <= (i - 1):
            n1 = elai.Element(i - 1, color)
            n2 = elai.Element(i, color)
            problem.setMatrix(n1, n2, 1)
        if (i + 1) < nnd:
            n1 = elai.Element(i, color)
            n2 = elai.Element(i + 1, color)
            problem.setMatrix(n1, n2, 1)
        if (i + 4) < nnd:
            n1 = elai.Element(i, color)
            n2 = elai.Element(i + 4, color)
            problem.setMatrix(n1, n2, 1)
        inhomogeneous.setVector(elai.Element(i, color), 1)

    A = problem.action()
    x = solution.ran()
    b = inhomogeneous.ran()
    direct(A, x, b)

if arg.afile is None or arg.bfile is None:
    demo()
    exit(0)

A = elai.ElaiMatrixDouble(elai.ProxyIfstream(arg.afile))
b = elai.ElaiVectorDouble(elai.ProxyIfstream(arg.bfile))
x = elai.ElaiVectorDouble(b.m())

imax = A.m() / 5
if arg.method == 'mumps':
    direct(A, x, b)
elif MPI.COMM_WORLD.Get_size() > 1:
    run_dist(A, x, b)
else:
    run(A, x, b)
