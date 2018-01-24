/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University.
Started: 9/28/2014 
Finished: N/A
*/

#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>

#include "cplex.h"

//FILE *inserted_data = NULL;


//struct branch_node
//{
//  	int seqnum;		// As given by CPLEX	
//	int split_dir;	 	// 0 for vertical, 1 for horizontal
//	double split_location;
//	double new_lb;		// For lower (or left) split
//	struct branch_node *prev;
//  	struct branch_node *next;
//};

struct branch_node
{
  	int seqnum;		// As given by CPLEX	
	int var_index;	 	// Index of variable to branch on
	double split_val;	// Previous value of variable (branches will use floor and ceiling of this value)
	struct branch_node *prev;
  	struct branch_node *next;
};

typedef struct branch_node branch_node;

struct pareto_branch_node
{
  	int seqnum;		// As given by CPLEX		 	
	double x_ub;
	double y_lb;
	double x_lb;
	double y_ub;
	struct pareto_branch_node *prev;
  	struct pareto_branch_node *next;
};

typedef struct pareto_branch_node pareto_branch_node;

struct mip_starts
{
  	int numstarts;		
  	int nzcnt;		 	
	int *beg;
	int *varindices;
	double *values;
	int *effortlevel;
};

typedef struct mip_starts mip_starts;

//struct user_data
//{
//	mip_starts *ws;
//	mip_starts *ob1;
//	mip_starts *ob2;
//	int ws_still_feas;
//	int ob1_still_feas;
//	int ob2_still_feas;
//	int ws_par_opt;
//	int ob1_par_opt;
//	int ob2_par_opt;
//	double ws_vals[2];
//	double ob1_vals[2];
//	double ob2_vals[2];
//	double *x_ws;
//	double *x1;
//	double *x2;
//};

//struct user_data
//{
//	double *x_ws;
//	double *x1;
//	double *x2;
//	int ws_still_feas;
//	int ob1_still_feas;
//	int ob2_still_feas;
//};

//typedef struct user_data user_data;

//void add_branch_node(int seqnum, int split_dir, double split_location, double new_lb);
void add_branch_node(int seqnum, int var_index, double split_val);

void add_pareto_branch(int seqnum, double x_ub, double y_lb, double x_lb, double y_ub);

void score_it(double *x);
void get_scores();

void free_coefs();
void free_xctype();

void pass_the_lps(CPXENVptr env, CPXLPptr prob1, CPXLPptr prob2);

void pass_a_mip(CPXLPptr prob1);

void close_all_env();

void PSA_all(CPXCENVptr env, CPXLPptr prob);

//extern double **stored_x;
//extern int x_rotation;
               
int CPXPUBLIC
nodeoperations (CPXCENVptr env,
                void       *cbdata,
                int        wherefrom,
                void       *cbhandle,
                int        *useraction_p);
                
int CPXPUBLIC 
usersetbranch (CPXCENVptr   env,
               void         *cbdata,
               int          wherefrom,
               void         *cbhandle,
               int          type,
               int          sos,
               int          nodecnt,
               int          bdcnt,
               const int    *nodebeg,
               const int    *indices,
               const char   *lu,
               const double *bd,
               const double *nodeest,
               int          *useraction_p);                
                
int CPXPUBLIC
userselectnode (CPXCENVptr env,
                void       *cbdata,
                int        wherefrom,
                void       *cbhandle,
                int        *nodenum_p,
                int        *useraction_p);
                
int CPXPUBLIC
userincumbent (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p);
               
int CPXPUBLIC
userincumbent_presolve (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p);
               
int CPXPUBLIC
userincumbent_cfe (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p);
               
int CPXPUBLIC
userincumbent2 (CPXCENVptr env3,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p);
               
//int CPXPUBLIC
//branch_nosol (CPXCENVptr env,
//           	void *cbdata,
//           	int wherefrom,
//          	void *cbhandle,
//           	int type,
//           	int sos,
//           	int nodecnt,
//           	int bdcnt,
//           	const int *nodebeg,
//           	const int *indices,
//           	const char *lu,
//           	const double *bd,
//           	const double *nodeest,
//           	int *useraction_p);
           	
//int CPXPUBLIC
//userheuristic (CPXCENVptr env,
//           void *cbdata,
//           int wherefrom,
//           void *cbhandle,
//           double *objval_p,
//           double *x,
//           int *checkfeas_p,
//           int *useraction_p);

int CPXPUBLIC
 cutcallback (CPXCENVptr env,
           void *cbdata,
           int wherefrom,
           void *cbhandle,
           int *useraction_p);
           
           int CPXPUBLIC 
usersetbranch3 (CPXCENVptr   env,
               void         *cbdata,
               int          wherefrom,
               void         *cbhandle,
               int          type,
               int          sos,
               int          nodecnt,
               int          bdcnt,
               const int    *nodebeg,
               const int    *indices,
               const char   *lu,
               const double *bd,
               const double *nodeest,
               int          *useraction_p);    
