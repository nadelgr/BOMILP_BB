The subdirectories contained herein hold the following files:

1) The original *.mps file as downloaded from MIPlib2017.
2) "original_instance.lp" -- A conversion of the above file to LP format. If the original objective was minimization, it has been converted to maximization.
3) "random_objective.lp" -- A new instance having all the same constraints as (2) but with a new randomly generated objective function.
4) "negative_objective.lp" -- A new instance having all the same constraints as (2) but whose objective function is the additive inverse of the objective in (2).

Note that in each *.lp file, all variables are given an objective coefficient, even if it is zero. This is to ensure that when these files are read into the BB procedure, the order in which variables are encountered is consistent for each file. BB *does not* function correctly if variables are encountered in different orders in different files. 
