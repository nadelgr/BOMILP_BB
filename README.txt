For instructions on compiling this software and to see examples of how to run the software from the terminal, see the INSTALL.txt file.

Notes:

(1) Setting Parameters: Before compiling, it may desireable for certain users to set various parameters in the file "user_set_parameters.h". That file also 	
	contains a description of each parameter and describes its function. By changing the name of the executable in the "Makefile," modifying these parameters,
	and re-compiling, one can create multiple executables, each with a different functionality.
    
(2) Output files: Running this software will produce two outputs:
	(a) The standard out displays human-readable updates throughout the execution of the software as well as a detailed summary of the results after termination.
	(b) A secondary output file is created, which can be named by the user (see INSTALL.txt), which contains the following information, tab delimited:
	
	1) Number of variables fixed by duality fixing
	2) Number of singleton columns discovered	
	3) Number of singleton columns fixed	
	4) Number of pairs of dominating columns discovered
	5) Number of bounds tightened using a dominated column
	6) Number of dominated columns fixed	
	7) Presolve Phase 1 time	
	8) Prepopulate time
	9) Presolve Phase 2 time
	10) Number of variable bounds tightened during ph2	
	11) Total BB time	
	12) Number of BB nodes explored
	13) Number of nodes at which local cuts were added
	14) Number of times branching was performed on an objective space disjunction
	15) Number of times a disjunctive cut was generated using an objective space disjucntion
	16) Number of MIPs solved during BB
	17) Number of BB nodes fathomed by PSA running to completion while maintaining integer-feasibility
	18) Number of BB nodes fathomed because the ideal point(s) of (subsets of) the dual bound was dominated
	19) Number of BB nodes fathomed because the ideal segment of the dual bound was dominated
	20) Number of BB nodes fathomed because the ideal point of a subset of the dual bound was dominated and the ideal segment of the remainder of the dual bound
	    was also dominated
	21) Number of BB nodes fathomed after each segment of the dual bound was generated, and each was dominated
	22) Number of variable bounds tightened during entire BB using probing before branching 
	23) Number of variable bounds tightened during probing due to infeasibility
	24) Number of variable bounds tightened during probing due to a dominated ideal point
	25) Number of variable bounds tightened during probing due to a dominated ideal segment
	26) Number of variable bounds tightened during probing due to a dominated lower bound set
	27) Total time spent in node processing function during BB
	28) Total time spent in user selected branching function during BB
	29) Total time spent probing during BB
	30) Total time spent generating local cuts during BB
	31) Total time spent generating disjcuntive cuts during BB 
	32) Total time spent solving MIPs during BB
	33) Maximum time spent running BB on a subset of the objective space when objective space splitting is employed
	34) Item (33) + Preprocessing/Presolve time (this should be approximately the time to run BB in parallel over each subset when exploiting objective space
	    gaps)
	35) Binary indicator of whether or not BB completely generated the Pareto set
	36) Total time spent until the start of duality gap calculation
	37) Number of open nodes at start of duality gap calculation
	38) Duality gap
	39) Duality gap - as a percent
	40) Length of nondominated subset of the dual bound
	41) Length of nondominated subset of the dual bound - measured as a percentage of the max of the f1-range and f2-range in the objective space
	
	This file is designed to have its contents copied and pasted into a spreadsheet so that multiple instances can be compared.
