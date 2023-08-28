## below is the name of the compiler
MCC     = mpicxx         #  mpiwrappper for the PGI C++ Compiler pgCC
##CFLAGS  = -fastsse -mcmodel=medium -tp k8-64e


CXXFLAGS = -fastsse -O3 -Mipa=fast -Minfo=all -mcmodel=medium
##CXXFLAGS = amd64a -fastsse -tp amd64 -mcmodel=medium 
OBJS    = EFIT.o

## name of the "target" file, or executable file of the compiling is 'test'

EFITsim: ${OBJS}  ## test is the name of the target/executable
	${MCC} ${CFLAGS} -o $@ ${OBJS}
##note that $@ means 'file name of the target'

EFIT.o:  EFIT.cpp   ## play.cpp is the name of the c++ file you're compiling
	${MCC} ${CFLAGS} -c $?

clean:
	-rm -f*.o *.bin core EFITsim a.out

##
##
## a.out is the default name of the executable file
##Makefiles explained
##test:		${OBJS} 
##	${MCC} ${CFLAGS} -o $@ ${OBJS} 
##  so ${OBJS}, ${MCC}, ${CFLAGS} are just variables names
##  what this line really looks like is
## test:	play.o
## 	 mpicxx  -fastsse -mcmodel=medium -tp k8-64e -o test play.o
##  so we are creating the executable 'test' by linking the object file 'play.o'
##  the object file or files as the case maybe are the 'dependent' files for 'test'
## play.o:		play.cpp
##	mpicxx  -fastsse -mcmodel=medium -tp k8-64e -c play.cpp
##
## here we are compiling an object file play.o  that is what the '-c' does 
## it tell the c++ compiler to make just the object file from the source code play.cpp
##   


