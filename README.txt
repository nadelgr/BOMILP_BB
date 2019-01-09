For instructions on compiling this software and to see examples of how to run the software from the terminal, see the INSTALL.txt file.

Notes:

(1) Setting Parameters: A variety of parameters can be modified through use of command line flags. See the end of this file for a list of such flags and their 
                        possible values.
    
(2) Output files: Running this software will produce two outputs:
	(a) The standard out displays human-readable updates throughout the execution of the software as well as a detailed summary of the results after termination.
	(b) If desired, a secondary output file is created (which can be named by the user using the appropriate command line flags). This file contains the 
	    following information, tab delimited:
	
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
	
	This file is written in "append" format, i.e., its contents are not overwritten each time a new instance is run. Thus, the file is designed to have its contents 
	copied and pasted into a spreadsheet so that multiple instances can be compared.
	
	_______________________________________________________________________________________________________________________________________________
	
	Command Line Flags: 
	
    1) -presolve1       Values: T - Presolve phase 1 is turned on
                                F - Presolve phase 1 is turned off (default)
                                
    2) -preprocessing   Values: T - Preprocessing is turned on (default) 
                                F - Preprocessing is turned off 
                                
    3) -ppType          Values: 0 - Automatically detect which preprocessing technique to use
                                1 - Use a preprocessing technique based on the weighted sum approach 
                                2 - Use a preprocessing technique based on the epsilon constraint approach (default)
                                3 - Use a hybrid preprocessing technique
                                
    4) -ppParam         Values: nonnegative integers -  This will adjust the number of solutions generated during preprocessing.
				                                        Increasing the value "should" increase the number of solutions generated. 
				                                        Hence, preprocessing will take longer, but could potentially improve 
				                                        fathoming throughout BB. Note that this value is coded into BOTH the 
				                                        epsilon-constraint and weighted-sum approaches and will therefore impact 
				                                        whichever procedure is indicated above. (default value - 1)
				                                        
	5) -ppExtra         Values: T - Utilize CPLEX's populate command instead of its mipopt command. This can help find extra 
	                                solutions during preprocessing for difficult problems, but it can also add significant time.
                                F - Use mipopt (default)
                                
    6) -ppCuts          Values: T - Each time a weighted sum problem is solved during preprocessing a local cut will be added along 
                                    the level curve associated with the optimal solution or best found dual value. This should serve 
                                    to improve the dual bound at the start of BB. However, if numerical difficulties arise during 
                                    the process of generating these cuts, this can result in integer feasible solutions being lost.
                                F - Do not add local cuts (default)
                                
    7) -ppGenetic       Values: nonnegative integers -  Indicates the number of times a simple genetic algorithm will be called 
							                            after each iteration of preprocessing in order to try to determine additional 
							                            nondominated integer feasible solutions not discovered during preprocessing. 
							                            (default value - 0)
							                            
	8) -presolve2       Values: T - Presolve phase 2 is turned on. Note that presolve phase 2 cannot be turned on if preprocessing 
	                                is turned off.
                                F - Presolve phase 2 is turned off (default)
                                
    9) -osFathom        Values: T - Turn on objective space fathoming (a.k.a. Pareto branching). (default)
                                F - Turn off objective space fathoming. 
                                
    10) -cutStrength    Values: 0 - Turn off CPLEX global cut generation.
                                1 - Use CPLEX default cut generation (default)
                                2 - Use CPLEX aggressive cut generation
                                
    11) -cutMulti       Values: T - Generate global cuts from a variety of weighted sum objectives prior to beginning BB. The goal is 
                                    to begin with the tightest dual bound possible.
                                F - Only generate global cuts from a single, pre-selected weighted-sum objective. (default)
                                
    12) -cutNum         Values: positive integers - If -cultMulti = T, the value of this parameter will indicate the number of 
                                                    iterations of global cut generation that will be attempted. Increasing the value 
                                                    also increases the number of weighted sum objectives used. (default value - 4)
    13) -cutLocal       Values: T - Attempt to generate local cuts at each node.
                                F - Do not attempt to generate local cuts at each node. (default)
    
    14) -intObj         Values: T - Attempt to exploit problem structure when one objective function depends only on integer variables. 
                                    The reason you may not want to do this is that it can be time consuming for some problems. For 
                                    example, if the coefficients on the integer variables are very small, although the Pareto solutions 
                                    are all found on a discrete number of verical (or horizontal) lines in the objective space, these 
                                    lines will be tightly packed and in this case working with only singletons in the objective space 
                                    may be too time consuming.
                                F - Do not Attempt to exploit problem structure in this way. (default)
    
    15) -bndRed         Values: T - Tighten variable bounds after processing each node, prior to branching. (default)
                                F - Do not tighten variable bounds after processing each node, prior to branching.
    
    16) -bndInf         Values: T - Limit the tightening attempt mentioned above to situations in which fixing a variable results in an 
                                    infeasible subproblem.
                                F - Do not limit the tightening attempt mentioned above. (default)
    
    17) -bndFull        Values: T - Generate the entire dual bound associated with a variable fixing when attempting to tighten bounds
                                    at each node.
                                F - Only tighten a bound when ideal points or segments are dominated. Bounds aren't tightened as often 
                                    this way, but generating the entire dual bound is typically quite time consuming. (default)
    
    18) -bndLimit       Values: T - Limit the number of variables for which bound tightening is attempted after processing each node. 
                                    The limitation is set to check 10% of the variables, beginning with a randomly selected index.
                                F - Do not limit the number of variables for which bound tightening is attempted. (default)
    
    19) -psaEarly       Values: T - Each time a segment forming the dual bound is generated, it is first extended to a segment that 
                                    dominates the remainder of the dual bound. If this extension is dominated, the node can be fathomed 
                                    without the need for generating any additional segments of the dual bound.
                                F - Each segment of the dual bound is generated until one is found that is non-dominated. (default)
    
    20) -psaEarlyN      Values: T - The method described above is employed, but only after "-psaIter" (below) iterations of the PSA.
                                F - The method described above is never employed. (default)
    
    21) -psaIter        Values: positive integers - See above. (default value - 15)
    
    22) -psaStop        Values: T - Stop generating line segments from the dual bound of a set of fixed integer variables when a line 
                                    segment is encountered that is dominated by the primal bound.
                                F - The entire dual bound of a set of fixed integer variables will be computed each time an integer 
                                    feasible solution is discovered. (default)
    
    23) -extraMIPs      Values: T - By default, during the cut callback we solve one single objective MIP and add a local cut along the 
                                    level curve in the objective space associated with the best known dual bound at the termination of 
                                    this MIP. This helps in strengthening the biobjective dual bound at the current node. Setting this 
                                    flag to 'T' causes three single objective MIPs to be solved instead of one.
                                F - Solve only one single objective MIP during the cut callback. (default)
    
    24) -gaps           Values: T - The solutions stored after preprocessing are scanned for large gaps between adjacent solutions. When 
                                    these gaps exist, the objective space is divided into subregions based on these gaps and a sub-BOMILP 
                                    is solved over each subregion. Note that splitting in this way often causes more nodes to be explored. 
                                    Nevertheless, computation time can be significantly reduced for some problems.
                                F - Traditional BB is conducted on the entire objective space at once. (default)
    
    25) -osDisj         Values: T - Each time a portion of the objective space is determined to be dominated which, if removed, results 
                                    in two disjoint remaining regions of the objective space, a polyhedral disjunction is created based on 
                                    these two regions and a valid cutting plane is generated for the convex hull of the union of the 
                                    disjoint sets. This cut is added to the model as a local cut. Note that this can be performed with or 
                                    without objective space fathoming.
                                F - No such cuts are generated. (default)
    
    26) -showProgress   Values: T - BB progress is printed to the screen in a similar fashion to CPLEX's output for single objective MIPs. 
                                    Each time "-showFreq" (below) nodes are processed, the following information is displayed: (i) the 
                                    number of nodes processed so far, (ii) the number of remaining nodes, and (iii) the current dual gap 
                                    (measured using the hausdorff metric - reported as a percentage of the maximum of the f_1 and f_2 
                                    ranges of the instance, or the hypervolume - also reported as a relative percentage, depending on the 
                                    value of the "-hausHyper" parameter below). Note that unlike single objective MIP, dual information is 
                                    not readily available and must be constructed in order to be reported, and so showing progress can 
                                    cause a significant increase in overall solution time.
                                F - Progress is not reported until BB terminates. (default)
    
    27) -ShowFreq       Values: positive integers - The number of nodes that should be processed in between each progress report. Note that
                                                    this value has no effect unless -showProgress = T. (default value - 5)
    
    28) -hausHyper      Values: A - Duality gaps are measured using the hausdorff metric.
                                Y - Duality gaps are measured using hypervolume. (default)
    
    29) -nodeSel        Values: T - CPLEX's default node selection is occasionally overridden. The default node selection in CPLEX is to 
                                    choose the node with best weighted-sum objective value as the next node to process. If this parameter 
                                    is set to 1, each time a new node needs to be selected, a random number is generated (uniform 0-1). If 
                                    the number is less than "-nodeSelPerc" (below), then 50% of the time the node with maximum f_1 value
                                    in its weighted-sum objective is selected for processing, and 50% of the time the node with maximum f_2 
                                    value in its weighted-sum objective is selected for processing. This can often allow more progress to 
                                    made in dual bound gap reduction throughout BB, although, when exhaustive BB is desired, overall number 
                                    of nodes processed and solution time does not seem to be significantly reduced, in general.
                                F - CPLEX default node selection is always used. (default)
    
    30) -nodeSelPerc    Values: reals in (0,1) - See above. (default value - 0.95)
    
    31) -nodeSelVar     Values: T - Instead of using "-nodeSelPerc" as a cut-off for whether or not the default CPLEX node selection scheme 
                                    is used, a value between "-nodeSelPerc1" and "-nodeSelPerc2" (both below) is used. Essentially, in this 
                                    case the value of "-nodeSelPerc" is overridden, initially set to "-nodeSelPerc1" and then increased by 
                                    0.005 each time a new node is selected, until the value reaches "-nodeSelPerc2". We allow for this 
                                    scheme because it seems that low (near 0) values of "-nodeSelPerc" result in more progress closing the 
                                    dual bound gap early in the BB procedure, while large (near 1) values result in more progress near the 
                                    end of BB. However, when exhaustive BB is desired, there does not seem to be a large impact on number 
                                    of nodes processed or overall solution time. Note that this parameter has no impact if -nodeSel = F.
                                F - Use "-nodeDelPerc" as described above. (default)
    
    32) -nodeSelPerc1   Values: reals in (0,1) - See above. (default value - 0.1)
    
    33) -nodeSelPerc2   Values: reals in (0,1) - See above. (default value - 0.95)
    
    35) -domCol         Values: T - The check for dominated columns during presolve is turned on.
                                F - The check for dominated columns during presolve is turned off. (default)
    
    36) -domColCuts     Values: T - If -domCol = T, the dominance relationship is used to generate disjunctive cuts while processing the 
                                    root node. For further information on this dominance relationship and the resulting disjunctions, see 
                                    "Progress in presolving for mixed integer programming," (Gamrath et al. 2015). 
                                F - No such cuts are generated. (default)
    
    37) -dualFix        Values: T - The check for duality fixing during presolve is turned on.
                                F - The check for duality fixing during presolve is turned off. (default)
    
    38) -singCol        Values: T - The check for singleton columns during presolve is turned on.
                                F - The check for singleton columns during presolve is turned off. (default)
    
    39) -limitType      Values: T - Single objective MIPs will be solved for a time up to the value of the parameter -mipTime. (default)
                                N - Single objective MIPs will be solved until -mipNode nodes have been processed. 
    
    40) -proveInfeas    Values: T - Single objective MIPs which return the status "time limit infeasible" will continue to be solved until 
                                    this status changes, even if the time -mipTime is exceeded. The motivation for enabling this feature is 
                                    that if the MIP being considered is truly infeasible, the MIPs at every child node will also be 
                                    infeasible. Thus, it is likely that we recieve the "time limit infeasible" at every child node until the 
                                    associated subproblem becomes easy enough to show that the MIP is infeasible within "time_limit." Hence,
				                    allowing the initial MIP to run until infeasibility is proved can significantly reduce overall BB time.
                                F - Single objective MIPs will cease to be processed after -mipTime seconds.  (default)
    
    41) -bbTime         Values: positive reals -    The maximum time to spend in BB (in seconds). Note that presolve and preprocessing will 
                                                    run to termination even if this limit is exceeded. (default value - 43200)
    
    42) -bbNodes        Values: positive integers - The maximum number of nodes to process during BB. (default value - 100000000)
    
    43) -mipTime        Values: positive reals -    The maximum time to spend on each single objective MIP solve during BB (in seconds). 
                                                    (default value - 60)
    
    44) -bbNodes        Values: positive integers - The maximum number of nodes to process for each single objective MIP solved during BB. 
                                                    (default value - 1000000000)
    
    45) -presolve1time  Values: positive reals -    The maximum amount of time to spend in presolve phase 1 (in seconds). 
                                                    (default value - 300)
    
    46) -presolve2time  Values: positive reals -    The maximum amount of time to spend in presolve phase 2 (in seconds). 
                                                    (default value - 20)
    
    47) -bndRedTime     Values: positive reals -    The maximum amount of time to spend tightening variable bounds after processing each 
                                                    node (in seconds). (default value - 1)
    
    48) -ppTime         Values: positive reals -    The maximum amount of time to spend in preprocessing (in seconds). 
                                                    (default value - 1800)
    
    49) -dualGapLimit   Values: positive reals -    BB is terminated when the dual gap falls below its value. Note that unlike single 
                                                    objective MIP, dual information is not readily available and must be constructed in 
                                                    order to be utilized. This information is NOT generated unless -showProgress = T. 
                                                    Therefore, unless -showProgress = T, this value has NO EFFECT. (default value - 0.3)
    
    50) -dualTime       Values: positive reals -    The maximum time to be spent attempting to generate the problem's dual bound 
                                                    (in seconds) if the overall time limit or node limit is exceeded. After this time limit 
                                                    is exceeded, the dual bound will be approximated. If twice this time limit is exceeded, 
                                                    BB will simply terminate, and indicate that the quality of the solution is unknown. 
                                                    (default value - 7200)
    
    51) -probePerc      Values: reals in [0,100] -  The percentage of integer variables which will be checked for probing. Probing has been 
                                                    shown to be effective, but can be extremely time consuming if there are a large number 
                                                    of integer variables. (default value - 100)
    
    52) -out            Values: a filename -    The name of the output file containing summary information. Note that this file is opened 
                                                in "append" format so that previous results are not overwritten. 
                                                (default value - bb_results.txt)
                                                Note: To suppress output files entirely, enter "none"
    
    53) -print          Values: T -             The Pareto set is displayed in the standard out after termination of BB.
                                F -             Output of the Pareto set is suppressed. (default)
                                a filename -    The Pareto set is printed in the specified file after termination of BB.
    
    54) -matlab         Values: T - If displayed, the Pareto set is pre-formatted so that it can be copied and pasted into the MATLAB file 
                                    "display_Pareto.m" (included), which creates a plot of these solutions. (default)
                                F - The Pareto set is displayed in the format: (i) Segments: (x1,y1) -- (x2,y2), (ii) Points: (x,y).
                                
