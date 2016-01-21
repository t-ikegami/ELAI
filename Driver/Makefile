CC      = mpic++
LD      = mpic++
#LD      = mpif90
#CFLAGS  = -W -O0 -std=c++11 -DELAI_USE_MUMPS -DELAI_USE_MPI -DELAI_PROFILE -DELAI_DEBUG
#CFLAGS  = -W -O3 -std=c++11 -DELAI_USE_MUMPS -DELAI_USE_MPI
CFLAGS  = -W -O0 -std=c++11 -DELAI_USE_MPI -DELAI_PROFILE -DELAI_DEBUG
LDFLAGS = 
#LDFLAGS = -nofor-main
INCLUDE = -I$(ELAI_INC)
#LIBS    = -L$(ELAI_LIB) -Wl,-Bstatic -ldmumps -lmumps_common -lscalapack -lmetis -lpord -Wl,-Bdynamic -lmkl_rt -limf -lifcore -lmpi_f90 -lmpi_f77 -lmpi_cxx
#LIBS    = -L$(ELAI_LIB) -Wl,-Bstatic -ldmumps -lmumps_common -lscalapack -lmetis -lpord -Wl,-Bdynamic -lmkl_rt -limf -lifcore -lmpi_mpifh -lmpi_usempi -lmpi_cxx
#LIBS    = -L$(ELAI_LIB) -Wl,-Bstatic -ldmumps -lmumps_common -lscalapack -lmetis -lpord -Wl,-Bdynamic -lmkl_rt -limf -lstdc++
LIBS    = -L$(ELAI_LIB) #-lstdc++
TARGET  = solv
OBJS    = solver.o

all: $(TARGET)
$(TARGET): $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

.cc.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

.PHONY: clean test1 test2 test3 demo
clean:
	rm -f $(TARGET) $(OBJS)
test1: $(TARGET)
	CTHRES=1e-08 FLEVEL=1 FTHRES=1e-08 STHRES=1e-05 KSP="SPBCGS" \
	mpiexec -n 1 $(TARGET) A1.mtx b1.mtx 2>&1
test2: $(TARGET)
	CTHRES=1e-08 FLEVEL=1 FTHRES=1e-10 STHRES=1e-05 KSP="SPBCGSA" \
	mpiexec -n 2 $(TARGET) A2.mtx b2.mtx 2>&1
test3: $(TARGET)
	CTHRES=1e-08 FLEVEL=1 FTHRES=1e-08 STHRES=1e-05 KSP="LU" \
	mpiexec -n 2 $(TARGET) A2.mtx b2.mtx 2>&1
test4: $(TARGET)
	CTHRES=1e-08 FLEVEL=1 FTHRES=1e-08 STHRES=1e-05 KSP="SPGMRES" \
	mpiexec -n 2 $(TARGET) A1.mtx b1.mtx 2>&1
demo: $(TARGET)
	CTHRES=1e-08 FLEVEL=1 FTHRES=1e-08 STHRES=1e-05 KSP="LU" \
	mpiexec -n 2 $(TARGET)
