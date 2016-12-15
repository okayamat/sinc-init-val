# Makefile ( n source file 1 execute file version )

PACKAGE	= SE_coll_example1
SRCS	= $(PACKAGE).cpp SE_trans.cpp
HEADS	= SE_trans.h
OBJS	= $(SRCS:.cpp=.o)
FILES	= SE_coll_example1 $(HEADS) $(SRCS)
VER	= `date +%Y%m%d`


### command and flags ###
# uncomment when debugging
#DEBUG	= -ggdb -pg # -lefence

# common (*.o)
LD	= g++-4.0
#CPPL_LD	= -L/sw/lib -llapack -lf77blas -latlas -lg2c -lm
CPPL_LD = -L/sw/lib -L/sw/lib/lapack/3.4.2 -lreflapack -lrefblas -lm
GSL_LD	= `/sw/bin/gsl-config --libs`
#-L/usr/lib/blas/atlas -lcblas  -L/usr/lib/lapack/atlas -llapack 
LDFLAGS	= -g -Wall -Wno-unknown-pragmas $(DEBUG) $(CPPL_LD) $(GSL_LD)
#LDLIBS	= -lm

# C (*.c)
CC	= gcc
OPTIMIZE= -O2
CFLAGS	= -g $(OPTIMIZE) -Wall $(DEBUG)
CPPFLAGS= -I.

# C++ (*.cpp)
CXX	= g++-4.0
GSL_FL	= `/sw/bin/gsl-config --cflags`
CPPL_FL	= -I/sw/include -I/sw/include/cpplapack -Wno-unknown-pragmas
CXXFLAGS= -g $(OPTIMIZE) -Wall $(DEBUG) $(CPPL_FL) $(GSL_FL) 

# Fortran77 (*.f)
FC	= f77
FFLAGS	= -Wall $(DEBUG)

# Pascal (*.p)
PC	= pc
PFLAGS	= -Wall $(DEBUG)

# etc
SHELL	= /bin/sh
RM	= rm -f
PROF	= gprof


### rules ###

.SUFFIXES:
.SUFFIXES: .o .c .cpp .f .p

all: $(PACKAGE)

$(PACKAGE): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $@ $(LDLIBS)

$(OBJS): $(HEADS) $(PACKAGE).mak

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@
.cpp.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
.f.o:
	$(FC) $(FFLAGS) -c $< -o $@
.p.o:
	$(PC) $(PFLAGS) $(CPPFLAGS) -c $< -o $@


### useful commands ###

clean:
	$(RM) $(PACKAGE) $(OBJS)
	$(RM) core gmon.out *~ #*#

tar:
	@echo $(PACKAGE)-$(VER) > .package
	@$(RM) -r `cat .package`
	@mkdir `cat .package`
	@ln $(FILES) `cat .package`
	tar cvf - `cat .package` | gzip -9 > `cat .package`.tar.gz
	@$(RM) -r `cat .package` .package

zip:
	zip -9 $(PACKAGE)-$(VER).zip $(FILES)


prof: run
	$(PROF) $(PACKAGE) | less

run: all
#	./$(PACKAGE) < sample-data | less
	./$(PACKAGE)
