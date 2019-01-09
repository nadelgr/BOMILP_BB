#
# Set these two directories
#

HOMEBIN=${HOME}/Dropbox\ \(Edinboro\ University\)/Edinboro/Research/github_bb/BOMILP_BB
CPLEX_DIR=$(HOME)/ILOG/CPLEX_Studio126/cplex

CC=gcc

#
# This can be left unchanged
#

SOLSTRUCT=tree

LDFLAGS= -O3 -g -L$(CPLEX_DIR)/lib/x86-64_linux/static_pic -lcplex -lpthread -lm -fno-stack-protector
CFLAGS=  -O3 -g -I$(CPLEX_DIR)/include/ilcplex -fno-stack-protector

OBJ = heuristics.o minor_functions.o max_$(SOLSTRUCT).o callbacks.o presolve_preprocessing.o biobjective_bb.o 
HEADER = heuristics.h minor_functions.h user_set_parameters.h max_$(SOLSTRUCT).h callbacks.h presolve_preprocessing.h bb-bicriteria.h biobjective_bb.h 

all: ${HOMEBIN}/bb_solver

clean:
	@rm *.o
	@[ -f ${HOMEBIN}/bb_solver ] 

${HOMEBIN}/bb_solver: $(OBJ) Makefile $(HEADER)
	@echo Linking $(@F)
	@$(CC) -DSOL_$(SOLSTRUCT) -o ${HOMEBIN}/bb_solver $(OBJ) $(LDFLAGS)

%.o: %.c Makefile $(HEADER)
	@echo [${CC}] $<
	@$(CC) -DSOL_$(SOLSTRUCT) $(CFLAGS) -c $< 
