#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>


//  This file is created for setting user controlled parameters. I will try my best to label each with what it does and its attainable values.
//  Some will be self explanatory. 

// ****************************** Parameters affecting performance ***********************************************

int presolve_phase_1 = 0; // Set to 1 to turn on presolve phase 1. Otherwise, set to 0.

int preprocessing = 1; // Set to 1 to turn on preprocessing. Otherwise, set to 0.

int preprocessing_method = 2; // Set to 0 in order to automatically detect which preprocessing technique to use, otherwise, set to 1 in order to use a
			      // preprocessing technique based on the weighted sum approach, set to 2 in order to use an approach based on the epsilon-contraint
			      // approach, or set to 3 to use a hybrid method. 
			      
int preprocessing_parameter = 1; // This parameter should be set to a non-negative integer. This will adjust the number of solutions generated during preprocessing.
				 // Increasing the value "should" increase the number of solutions generated. Hence, preprocessing will take longer, but could 
				 // potentially improve fathoming throughout BB. Note that this value is coded into BOTH the epsilon-constraint and weighted-sum
				 // approaches and will therefore impact whichever procedure is indicated above.

int search_for_extra_solutions_during_preprocessing = 0; // Set to 1 in order to utilize CPLEX's populate command instead of its mipopt command. This can help find
							 // extra solutions during preprocessing for difficult problems, but it can also add significant time.

int add_local_cuts_during_preprocessing = 0; // If this parameter is set to 1, each time a weighted sum problem is solved during preprocessing a local cut will be 
					     // added along the level curve associated with the optimal solution or best found dual value. This should serve to 
					     // improve the dual bound at the start of BB. If numerical difficulties arise during the process of generating these
					     // cuts, this can result in integer feasible solutions being lost.
					     
int times_to_try_heuristic_per_preprocessing_iteration = 0; // The value of this parameter indicates the number of times a simple genetic algorithm will be called 
							     // after each iteration of preprocessing in order to try to determine additional nondominated integer
							     // feasible solutions not discovered during preprocessing. Setting the parameter to 0 ensures that the
							     // genetic algorithm is never used. 

int presolve_phase_2 = 0; // Set to 1 to turn on presolve phase 2. Otherwise, set to 0. Note that presolve phase 2 is not employed unless preprocessing is used.

int objective_space_fathoming = 1; // Set to 1 if you desire objective space fathoming (a.k.a. Pareto branching). Otherwise, set to 0.

int CPLEX_cuts = 1; // Set to 1 if you want CPLEX to generate global cuts at default value. Set to 2 to have CPLEX generate cuts aggressively. To turn off CPLEX
		    // global cuts, set to 0.
		    
int generate_CPLEX_global_cuts_for_several_weighted_sum_objectives = 0; // Set to 1 to generate global cuts from a variety of weighted sum objectives prior to
									// beginning BB. The goal is to begin with the tightest dual bound possible. If the
									// parameter is set to 0, CPLEX will only ever generate global cuts from a single, 
									// pre-selected weighted-sum objective.

int num_iterations_weights = 4; // If the previous parameter is set to 1, the value of this parameter will indicate the number of iterations of global cut
				// generation that will be attempted. Increasing the value also increases the number of weighted sum objectives used.

int generate_local_cuts = 0; // Set to 1 if you want to spend time generating local cuts at each node. Otherwise, set to 0.

int integer_bb = 0; // Set to 1 if you want to exploit problem structure when one objective function depends only on integer variables. Otherwise, set to 0.
		    // The reason you may not want to do this is that it can be time consuming for some problems. For example, if the coefficients on the integer
		    // variables are very small, although the Pareto solutions are all found on a discrete number of verical (or horizontal) lines in the objective 
		    // space, these lines will be tightly packed and in this case working with only singletons in the objective space may be too time consuming.

int use_hausdorff = 0; // Set to 0 to perform exhaustive BB. Set to 1 to consider an element of the dual bound dominated by an element of the primal bound whenever
		       // the proximal hausdorff distance between the two is less than epsilon*max(SE_x - NW_x, NW_y - SE_y). If setting to 1, also choose a value 
		       // for epsilon.
double epsilon = 0.0001; // See above.

int bd_reduction_during_branching = 1; // Set to 1 if you want to tighten variable bounds after processing each node, prior to branching. Otherwise, set to 0.
int only_infeasibility = 0; // Set to 1 if you want to limit the tightening attempt to situations in which fixing a variable results in an infeasible subproblem.
			    // Otherwise, set to 0.
int check_bound = 0; // Set to 1 if you want to generate the entire dual bound associated with a variable fixing when attempting to tighten bounds. Otherwise,
    		     // set to 0. Setting to 0 will only tighten a bound when ideal points or segments are dominated. Bounds aren't tightened as often this way,
    		     // but generating the entire dual bound is typically quite time consuming.
int limit_bd_reduction = 0; // Set to 1 if you want to limit the number of variables for which bound tightening is attempted after processing each node. The 
			    // limitation is set to check 10% of the variables, beginning with a randomly selected index. Set to 0 if you want to attempt to 
			    // tighten the bounds of all integer variables. Note that regardless of the setting of this parameter, this process is still stopped
			    // after "max_time_bd_red_b4_branching" seconds. Set the value of "max_time_bd_red_b4_branching" below.
				     
int check_for_early_PSA_reduce = 0; // This parameter affects the performance of checking a dual bound for domination by the primal bound. If this parameter is 
				    // set to 1, each time a segment forming the dual bound is generated, it is first extended to a segment that dominates the 
				    // remainder of the dual bound. If this extension is dominated, the node can be fathomed without the need for generating any 
				    // additional segments of the dual bound. If the parameter is set to 0, each segment of the dual bound is generated until one is
				    // found that is non-dominated.
int check_for_early_PSA_reduce_after_n_iter = 0; // If this parameter is set to 1, the method described above is implemented after "n_iter" segments of the 
						 // dual bound have been generated. If setting to 1, also set "n_iter" below.
int n_iter = 15; // See above.

int stop_PSA_full_early_after_preprocessing = 0; // Set to 1 in order to stop generating line segments from the dual bound of a set of fixed integer variables when
						 // a line segment is encountered that is dominated by the primal bound. If set to 0, the entire dual bound of a
						 // set of fixed integer variables will be computed each time an integer feasible solution is discovered.

int solve_extra_mips_during_cut_callback = 0; // During the cut callback we solve one single objective MIP and add a local cut along the level curve in the 
					      // objective space associated with the best known dual bound at the termination of this MIP. This helps in 
					      // strengthening the biobjective dual bound at the current node. Setting this parameter to 0 leaves this process as
					      // is. Setting it to 1 causes three single objective MIPs to be solved instead of one.
					      
int exploit_objective_gaps = 0; // If this parameter is set to 1, the solutions stored after preprocessing are scanned for large gaps between adjacent solutions.
				// When these gaps exist, it may be beneficial to break the objective space into subregions based on these gaps and solve each
				// sub-BOMILP individually. Of course, splitting this way easily gives way for parallelization, which is clearly beneficial. 
				// Also note that splitting in this way often causes more nodes to be explored, though computation time can often be reduced.
				// If the parameter is set to 0, traditional BB is conducted on the entire objective space at once.
				
int generate_disjunctive_cuts_from_obj_space_disjunctions = 0; // If this parameter is set to 1, each time a portion of the objective space is determined to be
							       // dominated which, if removed, results in two disjoint remaining regions of the objective space, a
							       // polyhedral disjunction is created based on these two regions and a valid cutting plane is
							       // generated for the convex hull of the union of the disjoint sets. This cut is added to the model
							       // as a local cut. Note that this can be performed with or without objective space fathoming. If the
							       // parameter is set to 0 no such cuts are generated.

int show_progress = 0; // If this parameter is set to 1, BB progress is printed to the screen in a similar fashion to CPLEX's output for single objective MIPs. Each
		       // time "show_frequency" nodes (this parameter should be set below) are processed, the following information is displayed: (i) the number of 
		       // nodes processed so far, (ii) the number of remaining nodes, and (iii) the current dual gap (measured using the hausdorff metric -
		       // reported as a percentage of the maximum of the f_1 and f_2 ranges of the instance, or the hypervolume - also reported as a relative
		       // percentage, depending on the value of the "hausdorff_or_hypervolume" parameter). Note that unlike single objective MIP, dual
		       // information is not readily available and must be constructed in order to be reported, and so showing progress can cause a significant
		       // increase in overall solution time. If the parameter is set to 0, progress is not reported until BB terminates.

int show_frequency = 5; // If the parameter "show_progress" is set to 1, the value of this parameter determines the number of nodes that should be processed in
			 // between each progress report. If "show_progress" is set to 0, this parameter has no impact.
			 
int hausdorff_or_hypervolume = 2; // If this parameter is set to 1 duality gaps are measured using the hausdorff metric, if it is set to 2 duality gaps are
		                  // measured using the hypervolume measure.
			 
int control_node_selection = 0; // If this parameter is set to 1, CPLEX's default node selection is occasionally overridden. The default node selection in CPLEX is 
				// to choose the node with best weighted-sum objective value as the next node to process. If this parameter is set to 1, each time 
				// a new node needs to be selected, a random number is generated (uniform 0-1). If the number is below "control_percentage" (the 
				// value of this parameter should be set below) then 50% of the time the node with maximum f_1 value in its weighted-sum objective
				// is selected for processing, and 50% of the time the node with maximum f_2 value in its weighted-sum objective is selected for
				// processing. This can often allow more progress to made in dual bound gap reduction throughout BB, although, when exhaustive BB
				// is desired, overall number of nodes processed and solution time does not seem to be significantly reduced in general. If the
				// parameter is set to 0, CPLEX default node selection is always used. 
				
double control_percentage = .95; // See above description.

int allow_changing_control_percentage = 0; // If the parameter "control_node_selection" is set to 1, setting this parameter to 1 has the following impact. Instead
					   // of using "control_percentage" as a cut-off for whether or not the default CPLEX node selection scheme is used, a
					   // value between "starting_percentage" and "stopping_percentage" is used. Essentially, in this case the value of 
					   // "control_percentage" is overridden, initially set to "starting_percentage" and then increased by 0.005 each time a 
					   // new node is selected, until the value reaches "stopping_percentage." We allow for this scheme because it seems that 
					   // low (near 0) values of "control_percentage" result in more progress closing the dual bound gap early in the BB
					   // procedure, while large (near 1) values result in more progress near the end of BB. However, when exhaustive BB is 
					   // desired, there does not seem to be a large impact on number of nodes processed or overall solution time. Note that 
					   // this parameter has no impact is "control_node_selection" is set to 0.
					   
double starting_percentage = .1; // See above.

double stopping_percentage = .95; // See above.

int dominating_columns = 0; // If set to 0 the check for dominated columns during presolve will be turned off.

int generate_disjunctive_cuts_from_dominated_columns = 0; // If this parameter is set to 1, and dominating columns are discovered during presolve, the dominance 
							  // relationship is used to generate disjunctive cuts while processing the root node. For further 
							  // information on this dominance relationship and the resulting disjunctions, see "Progress in presolving
							  // for mixed integer programming," (Gamrath et al. 2015). If the parameter is set to 0, no such cuts are
							  // generated.

int duality_fixing = 0; // If set to 0 the check for duality fixing during presolve will be turned off.

int singleton_columns = 0; // If set to 0 the check for singleton columns during presolve will be turned off.

int time_vs_node_lim = 0; // If set to 0 single objective MIPs will be solved for a time up to the value of the parameter "time_limit." If set to 1, single 
			  // objective MIPs will be allowed to process at most "node_limit" nodes.
			  
int keep_solving_infeasible_MIPs = 0; // If set to 1, MIPs which return the status "time limit infeasible" will continue to be solved until this status changes,
				      // even if the time "time_limit" is exceeded. If set to 0, the MIP will cease to be processed after "time_limit" seconds. The
				      // motivation for ssetting this parameter to 1 is that if the MIP being considered is truly infeasible, the MIPs at every 
				      // child node will also be infeasible. Thus, it is likely that we recieve the "time limit infeasible" at every child node
				      // until the associated subproblem becomes easy enough to show that the MIP is infeasible within "time_limit." Hence,
				      // allowing the initial MIP to run until infeasibility is proved can significantly reduce overall BB time.

// ****************************** Limits  ******************************************************************

double max_time = 43200.; // Set the maximum time to spend in BB. Note that presolve and preprocessing will run to termination even if this limit is exceeded.

int max_nodes = 100000000; // Set the maximum number of nodes to process during BB. 

double time_limit = 60.; // Set the maximum time to spend on each single objective MIP solve during BB.

int node_limit = 1000000000; // Set the maximum number of nodes to process for each single objective MIP solved during BB.

double max_time_phase_2 = 20.; // Set the maximum amount of time to spend in presolve phase 2.

double max_time_bd_red_b4_branching = 1.; // Set the maximum amount of time to spend tightening variable bounds after processing each node.

double max_preprocessing_time = 1800.; // Set the maximum amount of time to spend in preprocessing.

double max_phase1_time = 300.; // Set the maximum amount of time to spend in presolve phase 1.

double duality_gap_limit = .3; // If this parameter is set to a non-zero value, BB is terminated when the dual gap falls below its value. Note that unlike single
			        // objective MIP, dual information is not readily available and must be constructed in order to be utilized. This information is
			        // NOT generated unless the value of the parameter "show_progress" is set to 1. Therefore, unless "show_progess" is set to 1, the
			        // value of this parameter has NO EFFECT.
			        
double max_time_build_dual_bound = 7200.; // Set the maximum time the code should spend attempting to generate the problem's dual bound if the time limit or node
					  // limit is exceeded. After this time limit is exceeded, the dual bound will be approximated. If 2x this time limit is
					  // exceeded, BB will simply terminate, and indicate that the quality of the solution in unknown.

double percentage_of_integer_variables_to_check_for_probing_during_presolve = 100.; // Set to a value between 0 and 100. This will indicate the percentage of
										    // integer variables which will be checked for probing. Probing has been shown
										    // to be effective, but can be extremely time consuming if there are a large
										    // number of integer variables.
