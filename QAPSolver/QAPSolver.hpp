/*****************************************************************/
// Implementation of the breakout local search of: 
// "Breakout local search for the quadratic assignment problem", 
// Applied Mathematics and Computation 219(9), 1991, 4800-4815.
//
// Data file format: 
//  n,
// (nxn) flow matrix,
// (nxn) distance matrix
//
// Copyright : Una Benlic, 2012
// This code can be freely used for non-commercial purpose.
// Any use of this implementation or a modification of the code
// must acknowledge the work of U. Benlic and J.K. Hao
/****************************************************************/

#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <map>
#include<ctime>
#include<cmath>

#include <vector>

namespace qap {
#ifdef _DEBUG
#define time_limit 30
#else
#define time_limit 60 * 10
#endif // _DEBUG

#define mem_size 100000
typedef int*   type_vector;
typedef long** type_matrix;

const long infinite = 999999999;

using namespace std;

long iteration = 0;

//BLS parameters
double r1 = 0.7,r2 = 0.2; // Parameters used for directed perturbation. Number of moves elapsed before perturb. 
double init_pert_str = 0.15; //initial perturbation strength
double perturb_str; // actual perturbation strength
int T = 2500;
double P0=0.75, Q=0.3;


void transpose(int & a, int & b) {long temp = a; a = b; b = temp;}


/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
long compute_delta(int n, type_matrix & a, type_matrix & b,
                   type_vector & p, int i, int j)
 {long d; int k;
  d = (a[i][i]-a[j][j])*(b[p[j]][p[j]]-b[p[i]][p[i]]) +
      (a[i][j]-a[j][i])*(b[p[j]][p[i]]-b[p[i]][p[j]]);
  for (k = 1; k <= n; k = k + 1) if (k!=i && k!=j)
    d = d + (a[k][i]-a[k][j])*(b[p[k]][p[j]]-b[p[k]][p[i]]) +
            (a[i][k]-a[j][k])*(b[p[j]][p[k]]-b[p[i]][p[k]]);
  return(d);
 }

/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
long compute_delta_part(type_matrix & a, type_matrix & b,
                        type_vector & p, type_matrix & delta, 
                        int i, int j, int r, int s)
  {return(delta[i][j]+(a[r][i]-a[r][j]+a[s][j]-a[s][i])*
     (b[p[s]][p[i]]-b[p[s]][p[j]]+b[p[r]][p[j]]-b[p[r]][p[i]])+
     (a[i][r]-a[j][r]+a[j][s]-a[i][s])*
     (b[p[i]][p[s]]-b[p[j]][p[s]]+b[p[j]][p[r]]-b[p[i]][p[r]]) );
  }
void update_matrix_of_move_cost(int i_retained, int j_retained,long n, type_matrix & delta, type_vector & p, type_matrix & a, type_matrix & b)
{
     int i, j;
     for (i = 1; i < n; i = i+1) for (j = i+1; j <= n; j = j+1)
        if (i != i_retained && i != j_retained && 
            j != i_retained && j != j_retained)
         {delta[i][j] = 
            compute_delta_part(a, b, p, delta, 
                               i, j, i_retained, j_retained);}
        else
         {delta[i][j] = compute_delta(n, a, b, p, i, j);};
}

void apply_move(type_vector & p,long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, 
		int i_retained, int j_retained, clock_t & end_time)
{
	if(i_retained!=-1 && j_retained!=-1) // apply the selected perturbation move
	{
		last_swaped[i_retained][j_retained] = iteration;
                last_swaped[j_retained][i_retained] = iteration;
		swap(p[i_retained], p[j_retained]);
		current_cost = current_cost + delta[i_retained][j_retained];
		update_matrix_of_move_cost(i_retained, j_retained, n, delta, p, a, b);
		if(current_cost<best_cost)
		{
			best_cost = current_cost;
                        end_time = clock();
			iter_without_improvement = 0;
			for (int m = 1; m <= n; m = m+1) best_sol[m] = p[m];
			//cout << "Solution of value " << best_cost<< " found. " <<iteration<< endl;

		}
	} iteration++;
}
void tabu_search_perturb(type_vector & p, long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, long init_cost, 
                clock_t & end_time)
{
	int i_retained=-1, j_retained=-1,  i, j, min_delta;
	min_delta = infinite;
	min_delta = 999999999;
        for (i = 1; i < n; i = i + 1) 
        {
           for (j = i+1; j <= n; j = j+1)
           {
             if((current_cost + delta[i][j])!=init_cost && delta[i][j] < min_delta && ((last_swaped[i][j] +n*r1+ rand()%(static_cast<int>(n*r2)))<iteration or (current_cost + delta[i][j])<best_cost))
             {
                i_retained = i; j_retained = j;
                min_delta = delta[i][j];
             };
           };
        };
	apply_move(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, end_time);  
}
void recency_based_perturb(type_vector & p, long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, clock_t & end_time)
{
	int i, j, i_retained, j_retained, min = infinite;

        for (i = 1; i < n; i = i + 1)
        { 
           for (j = i+1; j <= n; j = j+1)
           {
              if(last_swaped[i][j]<min)
              {
                 i_retained = i; j_retained = j;
                 min = last_swaped[i][j];
              };
           };
        };
	apply_move(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, end_time); 
}
void random_perturb(type_vector & p, long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, long init_cost, clock_t & end_time)
{
	int i_retained = 1+rand()%n;
        int j_retained = 1+rand()%n;
        if(i_retained>j_retained)
             swap(i_retained, j_retained);
        while(i_retained==j_retained || (current_cost + delta[i_retained][j_retained])==init_cost)
        {
             j_retained = 1+rand()%n;
             if(i_retained>j_retained)
               swap(i_retained, j_retained);
        }
	apply_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, end_time); 
}
void perturbe(type_vector & p,long n, type_matrix & delta, long & current_cost, type_matrix & a, type_matrix & b,
             type_matrix & last_swaped,  int & iter_without_improvement, type_vector & best_sol, long & best_cost, clock_t & end_time)
{
   int i_retained=-1, j_retained=-1, k, i, j, min_delta, min;
   long cost = current_cost; 

   double d = static_cast<double>(iter_without_improvement)/T;
   double e = pow(2.718, -d);
   int perturbation_bit = 0;
   if(e<P0)
      e=P0;
   if(e>(static_cast<double>(rand()%101)/100.0))
        perturbation_bit=1;
 
   for(k =0; k<(perturb_str); k++)
   {
      i_retained=-1; j_retained=-1;
      if(perturbation_bit==1)
      {
           tabu_search_perturb(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, cost, end_time);
      }
       else
       {  
         if(Q>(static_cast<double>(rand()%101)/100.0))
           recency_based_perturb(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, end_time);
         else
            random_perturb(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, cost, end_time);
          
       }
      
   }
}

void allocate_memory_and_initialize(type_vector & p, long n, type_matrix & delta, type_matrix & last_swaped, long & current_cost, 
		type_matrix a, type_matrix b, type_vector & best_sol, long & best_cost)
{
	int i, j;
	/***************** dynamic memory allocation *******************/
	p = new int[n+1];
	delta = new long* [n+1];
	last_swaped = new long* [n+1];
	for (i = 1; i <= n; i = i+1)
	{ 
		delta[i] = new long[n+1]; 
		last_swaped[i] = new long [n+1];
	}
	/************** current solution initialization ****************/
	for (i = 1; i <= n; i = i + 1) 
		p[i] = best_sol[i];
	/************** tabu list initialization ****************/
	for (i = 1; i <= n; i = i+1)
	{
		for(j=1; j<=n; j++)
		{
			last_swaped[i][j] = 0;
		}
	}
	/********** initialization of current solution value ***********/
	/**************** and matrix of cost of moves  *****************/
	current_cost = 0;
	for (i = 1; i <= n; i = i + 1) 
		for (j = 1; j <= n; j = j + 1)
		{
			current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
			if (i < j) 
			{
				delta[i][j] = compute_delta(n, a, b, p, i, j);
			};
		};
	best_cost = current_cost;

}

bool best_improvement_move(type_vector & p,long n, type_matrix & delta, long & current_cost, type_matrix a, type_matrix b,
		type_matrix & last_swaped, int & iter_without_improvement, long & best_cost, type_vector & best_sol, clock_t & end_time)
{
    int i, j, i_retained, j_retained;
    long min_delta = infinite;   // retained move cost
    //select the best swap move for the descent local search phase 
    for (i = 1; i < n; i = i + 1)
    { 
		for (j = i+1; j <= n; j = j+1)
		{
			if(delta[i][j] < min_delta)
			{
				i_retained = i; j_retained = j;
				min_delta = delta[i][j];
			};
		};
    };
    // apply the above selected move if it leads to a solution better than the current one
    if((current_cost + delta[i_retained][j_retained])<current_cost)
    {
	apply_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, i_retained, j_retained, end_time); 
		return 1;
    }else
		return 0;
}  
void determine_jump_magnitude(int & iter_without_improvement, int descent_num, long previous_cost,
		long current_cost, long n)
{
	//the following lines determine the number of perturbation moves (jump magnitude)
	if(iter_without_improvement>T) // if the best found solution is not improved during the last T descent phases 
	{
		// strong pertubation required
		iter_without_improvement = 0;
		perturb_str = n*(0.4 + (static_cast<double>(rand()%20)/100.0));
	}
	else if((descent_num!=0 && previous_cost != current_cost) || (descent_num!=0  && previous_cost != current_cost)) 
		// Escaped from the previous local optimum.
	{
		iter_without_improvement++;
		perturb_str=ceil(n*init_pert_str); if(perturb_str<5) perturb_str = 5;
	}
	else if(previous_cost == current_cost) // search returned to the previous local optimum
	{
		perturb_str+= 1;
	}
}     
void BLS(long n,                  // problem size
		type_matrix  a,         // flows matrix
		type_matrix  b,         // distance matrix
		type_vector & best_sol,  // best solution found
		long & best_cost,      // cost of best solution
                long num_iterations, long best_objective, clock_t & end_time)

{
	type_vector p;                        // current solution
	type_matrix delta;                   
	type_matrix last_swaped; // keeps track of the iteration number when a move was last performed

	long current_iteration;              
	long current_cost, previous_cost;    // current sol. value and previous sol. value
	int descent_num;
	int iter_without_improvement = 0; // counter of the number of consecutive descent phases with no improvement of the best solution 
	perturb_str = ceil(init_pert_str*n); // initialize the number of perturbation moves
	bool descent;

	allocate_memory_and_initialize(p, n, delta, last_swaped, current_cost, a, b, best_sol, best_cost);
	previous_cost = current_cost;
        clock_t start = clock();
	/******************** main BLS loop ********************/
	for (current_iteration = 1; current_iteration <= num_iterations; current_iteration ++)
	{
		//best improvement descent procedure
		descent = best_improvement_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_cost, best_sol, end_time);
                descent_num=0;
		while(descent) // while an improving solution in the neighborhood exists
                {
		   descent = best_improvement_move(p, n, delta, current_cost, a, b, last_swaped, iter_without_improvement, 
                                                   best_cost, best_sol, end_time);      
                   descent_num++;
                }

		// a this point, a local optimum is attained and BLS enters a perturbation phase
		determine_jump_magnitude(iter_without_improvement, descent_num, previous_cost, current_cost, n);

                //update the objective value of the previous local optimum
                previous_cost = current_cost;
        
                perturbe(p,n, delta, current_cost, a, b, last_swaped, iter_without_improvement, best_sol, best_cost, end_time);
                 if(ceil((clock() - start)/static_cast<double>(CLOCKS_PER_SEC))>=time_limit or best_cost==best_objective)
                   break;
	} 
        end_time-=start;
	// free memory
	delete[] p;
	for (int i=1; i <= n; i = i+1)
	{ 
		delete[] delta[i];
		delete[] last_swaped[i];
	}
	delete[] delta;  delete[] last_swaped; 
} // BLS

void generate_random_solution(long n, type_vector  & p)
{
  int i;
  for (i = 0; i <= n; i = i+1) p[i] = i;
  for (i = 1; i <  n; i = i+1) transpose(p[i], p[i + rand()%(n-i+1)]);
}

void load_problem( int &n, type_matrix &a, type_matrix &b, long & best_objective)
{
	cin >> best_objective >> n; 
	a = new long* [n+1];
	b = new long* [n+1];
	for (int i = 1; i <= n; i = i+1) 
	{
		a[i] = new long[n+1];
		b[i] = new long[n+1];
	}

	for (int i = 1; i <= n; i = i+1) 
		for (int j = 1; j <= n; j = j+1)
			cin >> a[i][j];
	for (int i = 1; i <= n; i = i+1) 
		for (int j = 1; j <= n; j = j+1)
			cin >> b[i][j];
}
void load_problem_from_datfile(int &n, type_matrix &a, type_matrix &b, long & best_objective)
{
	ifstream ifs("Instance/qapdata/tai20a.dat");
	if (!ifs.is_open()) { return; }
	
	best_objective = 703482;
	ifs >> n;
	a = new long*[n + 1];
	b = new long*[n + 1];
	for (int i = 1; i <= n; i = i+1)
	{
		a[i] = new long[n + 1];
		b[i] = new long[n + 1];
	}

	for (int i = 1; i <= n; i = i+1)
		for (int j = 1; j <= n; j = j+1)
			ifs >> a[i][j];
	for (int i = 1; i <= n; i = i+1)
		for (int j = 1; j <= n; j = j+1)
			ifs >> b[i][j];
}

int n;                    // problem size
type_matrix a, b;         // flows and distances matrices
type_vector solution;     // solution (permutation) 
long cost, best_objective;                // solution cost
int i, j; //variables for iterations

void run_bls(int &n, type_matrix &a, type_matrix &b, long &best_objective) {
	srand(time(NULL));
	/****************** dynamic memory allocation ******************/
	solution = new int[n + 1];
	/************** generate random solution and run BLS algorithm **************/

	long num_iterations = 2000000000;
	clock_t time;

	generate_random_solution(n, solution);
	iteration = 0;
	BLS(n, a, b, solution, cost, num_iterations, best_objective, time);

	//cout << "Solution cost " << cost << " " << (time) / static_cast<double>(CLOCKS_PER_SEC) << " " << 100.0*static_cast<double>(cost - best_objective) / static_cast<double>(best_objective) << endl;
	//cout << "Solution detail "; for (int i = 1; i <= n; i = i + 1) { cout << solution[i] << " "; } cout << endl;
}

void test_qap() {
    //load_problem(n, a, b, best_objective);
	load_problem_from_datfile(n, a, b, best_objective);
	
	run_bls(n, a, b, best_objective);
	
	delete[] solution;
	for (int i = 1; i <= n; i = i + 1) {
		delete[] a[i];
		delete[] b[i];
	}
	delete[] a; delete[] b;
}

void run_qap(const vector<vector<int>> &flow_matrix, const vector<vector<int>> &distance_matrix, vector<int> &sol) {
	// load_problem_from_vector
	best_objective = 0;
	n = flow_matrix.size();
	a = new long*[n + 1];
	b = new long*[n + 1];
	for (int i = 1; i <= n; i = i + 1) {
		a[i] = new long[n + 1];
		b[i] = new long[n + 1];
	}
	for (int i = 1; i <= n; i = i + 1)
		for (int j = 1; j <= n; j = j + 1)
			a[i][j] = flow_matrix[i - 1][j - 1];
	for (int i = 1; i <= n; i = i + 1)
		for (int j = 1; j <= n; j = j + 1)
			b[i][j] = distance_matrix[i - 1][j - 1];

	run_bls(n, a, b, best_objective);

	sol.clear(); sol.reserve(n);
	for (int i = 1; i <= n; i = i + 1) { sol.push_back(solution[i]); }

	delete[] solution;
	for (int i = 1; i <= n; i = i + 1) {
		delete[] a[i];
		delete[] b[i];
	}
	delete[] a; delete[] b;

	// check cost
	//int check_cost = 0;
	//for (int i = 0; i < n; ++i)
	//	for (int j = 0; j < n; ++j)
	//		check_cost += flow_matrix[i][j] * distance_matrix[sol[i] - 1][sol[j] - 1];
	//cout << (cost == check_cost) << endl;
}

}
