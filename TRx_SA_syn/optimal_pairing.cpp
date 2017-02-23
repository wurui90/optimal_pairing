//Simulated anealing algorithm

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <utility>
#include <random>
#include <assert.h>
#include <ctime>

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include "../shared/optimal_pairing_shared.cpp"
//#include "../shared/utility.cpp"

#define CHAN_NUM 5
#define PRINT_DB(X) printf("X: %lf\n", X);

using namespace std;


int main(int argc, char** argv)
{
    //Configuration parameters
    int TRX_NUM = stoi(argv[2]);
	int spacingG = stoi(argv[1]);
    double channel_spacing = 0.64*spacingG/80; //0.64 nm for 80 GHz channel spacing
	double lambda1 = 1.0;
//  double wl0 = 1310; //center wavelength
//  double wl_std = 0.2;
	
	//Read the WLs
	char WL_file_name[100];
	sprintf(WL_file_name, "../syn_WLs_%dG_num_%d.csv", spacingG, TRX_NUM);
	printf("Open file %s\n", WL_file_name);
	vector< vector<double> > TRx_WLs = ReadWLFromFile(WL_file_name, TRX_NUM);

	//Generate the cost matrix
	vector< vector<double> > costMat = GenerateCostMat(TRx_WLs, channel_spacing);

	//The greedy algorithm
	//Use its result as the initial assignment for the SA algorithm
	vector< pair<int, int>> assignment = AssignGreedy(costMat);
	PrintAssignDetail(assignment, costMat);
	double initCost = AssignCost(assignment, costMat, lambda1);

	//The simulated annealing algorithm
	int pair_num = assignment.size();
	double cost0 = initCost; //the previous-state cost
	double cost1 = 0; //the current-state cost
	
	double kT = initCost*0.5;		//Initial temperature for the SA algorithm
	double kT_final = initCost*0.001; //Final temperature
	double kT_scale_rate = 0.9;		//Cooling rate
	int loops_per_temp = 8000;		//Loop number for each temperature to reach equilibrium
	double E_min_per_temp0 = cost0;
	double E_min_per_temp1 = cost0;
	double E_min_global = cost0;
	vector< pair<int, int> > assignment_optimal(assignment);
	double E_change_criteria = 0.01;
	
    srand((int)time(0));
	int perturb_cnt = 0;
	while(kT >= kT_final)
	{
		E_min_per_temp1 = cost0;

		//loop for one temperature
		for(int i = 0; i < loops_per_temp; i++)
		{
			//Perturb the system
			perturb_cnt++;
			printf("%d:\t", perturb_cnt);
			//Randomly select two pairs to shuffle
			int pair1 = rand()%pair_num;
			int pair2 = rand()%pair_num;
			printf("kT: %lf \t", kT);
			printf("Shuffle pair %d and %d:  ", pair1, pair2);
			//Select shuffle pattern
			int shuffle_pattern = rand()%2;
			int tmp = 0;
			if(shuffle_pattern == 0)
			{
				tmp = assignment[pair1].first;
				assignment[pair1].first = assignment[pair2].first;
				assignment[pair2].first = tmp;

			}
			else
			{
				tmp = assignment[pair1].first;
				assignment[pair1].first = assignment[pair2].second;
				assignment[pair2].second = tmp;
			}
			
			//Compute delta_E
			cost1 = AssignCost(assignment, costMat, lambda1);
			double delta_E = cost1 - cost0;
			if(delta_E < 0)
			{
				//accept the perturbation
				E_min_per_temp1 = cost1;
				printf("Accept\t");
			}
			else
			{
				double rand01 = ((double)rand())/RAND_MAX;
				if(rand01 < exp(-delta_E/kT))
				{
					printf("Accept\t");
					//accept the perturbation
				}
				else
				{
					printf("Reject\t");
					//reject the perturbation 
					cost1 = cost0;
					if(shuffle_pattern == 0)
					{
						tmp = assignment[pair1].first;
						assignment[pair1].first = assignment[pair2].first;
						assignment[pair2].first = tmp;

					}
					else
					{
						tmp = assignment[pair1].first;
						assignment[pair1].first = assignment[pair2].second;
						assignment[pair2].second = tmp;
					}
		
				}
			}
			
			printf("current cost: %.2lf\tUptodate min cost:%.2lf\n", cost1, E_min_global);
			if(cost1 < E_min_global) {E_min_global = cost1; assignment_optimal = assignment;}
			cost0 = cost1;
		}
	
		//Determine whether we continue to cool down and exit
		//If the minimum cost or energy during the last two temperatures are still changing
	/*	if(abs(E_min_per_temp1 - E_min_per_temp0)/E_min_per_temp0 > E_change_criteria)
		{
			;
		}
		else break;
	*/		
		kT = kT*kT_scale_rate;
		E_min_per_temp0 = E_min_per_temp1;

	}

	printf("Initial cost: %lf\n", initCost);
	printf("Global min cost: %lf\n", E_min_global);
	printf("Optimal assignment:\n");
	vector<double> optimal_array = PrintAssignDetail(assignment_optimal, costMat);	

	//print the nearest pairing
	vector< pair<int, int>> nearest_pairing;
	for(int i = 0; i < floor(TRX_NUM/2); i++)
	{
		nearest_pairing.push_back(pair<int, int> (2*i, 2*i + 1));
	}
	printf("Nearest pairing:\n");
    vector<double> tuning_cost_array_nearest = PrintAssignDetail(nearest_pairing, costMat);
	printf("Mean: %f, std: %f\n", Vmean(tuning_cost_array_nearest), Vstd(tuning_cost_array_nearest));
	printf("power saving %f:\n", 1 - Vmean(optimal_array)/Vmean(tuning_cost_array_nearest));
	printf("std reduction %f:\n", 1 - Vstd(optimal_array)/Vstd(tuning_cost_array_nearest));
    return 0;	
}


    //Randomly generate the WLs
/*    normal_distribution<double> norm_dist(0.0, wl_std);
    default_random_engine rng;

    printf("Synthetic WLs:\n");

    for(int i = 0; i < TRX_NUM; i++)
    {
        printf("Device#%d:", i);
        for(int j = 0; j < CHAN_NUM; j++)
        {
            TRx_WLs[i][j] = wl0 + channel_spacing*j + norm_dist(rng);
            TRx_WLs[i][j + CHAN_NUM] = wl0 + channel_spacing*j + norm_dist(rng);
        
            printf("%f\t", TRx_WLs[i][j]);
        }
        printf("\n");
    }
*/


