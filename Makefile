#
# Set these two directories
#

HOMEBIN=${HOME}/Desktop/tree_version
CPLEX_DIR=$(HOME)/ILOG/CPLEX_Studio126/cplex

CC=gcc

#
# This can be left unchanged
#

SOLSTRUCT=tree

LDFLAGS= -O3 -L$(CPLEX_DIR)/lib/x86-64_linux/static_pic -lcplex -lpthread -lm
CFLAGS=  -O3 -I$(CPLEX_DIR)/include/ilcplex

OBJ = max_$(SOLSTRUCT).o callbacks.o presolve_preprocessing.o biobjective_bb.o 
HEADER = user_set_parameters.h max_$(SOLSTRUCT).h callbacks.h presolve_preprocessing.h bb-bicriteria.h biobjective_bb.h 

all: ${HOMEBIN}/bb_solver

clean:
	@rm *.o
	@[ -f ${HOMEBIN}/bb_solver ] 

${HOMEBIN}/bb_solver_time: $(OBJ) Makefile $(HEADER)
	@echo Linking $(@F)
	@$(CC) -DSOL_$(SOLSTRUCT) -o ${HOMEBIN}/bb_solver $(OBJ) $(LDFLAGS)

%.o: %.c Makefile $(HEADER)
	@echo [${CC}] $<
	@$(CC) -DSOL_$(SOLSTRUCT) $(CFLAGS) -c $< 
