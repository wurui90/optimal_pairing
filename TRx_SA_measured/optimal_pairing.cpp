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

#define CHAN_NUM 5
#define PRINT_DB(X) printf("X: %lf\n", X);

using namespace std;

//tuning distance of aligning wl1 and wl2 with given wavelength grid
double TuningDistWL(vector<double> wl1, vector<double> wl2, double channel_spacing)
{
    double left_WL = -1.0;
//    printf("%d\n", (int)wl1.size());
    assert(wl1.size() == CHAN_NUM);
    assert(wl2.size() == CHAN_NUM);
    for(int i = 0; i < wl1.size(); i++)
    {
        left_WL = max(left_WL, wl1[i] - channel_spacing*i);
        left_WL = max(left_WL, wl2[i] - channel_spacing*i);
    }
    
    double ret = 0;
    for(int i = 0; i < wl1.size(); i++)
    {
        assert(left_WL + i*channel_spacing >= wl1[i]);
        assert(left_WL + i*channel_spacing >= wl2[i]);
        ret += (left_WL + i*channel_spacing - wl1[i]) + (left_WL + i*channel_spacing - wl2[i]);
    }

    return ret;
}

//Tuning cost with TRx1 and TRx2, each with 5x2 = 10 wavelengths 
double TuningDistTRx(vector<double> wl1, vector<double> wl2, double channel_spacing)
{
   // return abs(wl1[0] - wl2[0]);
    return TuningDistWL(vector<double> (wl1.begin(), wl1.begin() + CHAN_NUM), vector<double> (wl2.begin() + CHAN_NUM, wl2.end()), channel_spacing) + TuningDistWL(vector<double> (wl1.begin() + CHAN_NUM, wl1.end()), vector<double> (wl2.begin(), wl2.begin() + CHAN_NUM), channel_spacing);
}

//Calculate the cost for the given assignment
double AssignCost(const vector< pair<int, int> >& assignment, const vector< vector<double> >& costMat, double lambda1)
{
	double A_sum = 0;
	double A_sqr_sum = 0;
	int pair_num = assignment.size();
	for(int i = 0; i < assignment.size(); i++)
	{
		A_sum += costMat[assignment[i].first][assignment[i].second];
		A_sqr_sum += costMat[assignment[i].first][assignment[i].second]*costMat[assignment[i].first][assignment[i].second];
	}
	return A_sum/pair_num + lambda1*sqrt(A_sqr_sum/pair_num - A_sum*A_sum/pair_num/pair_num);
	
}

void PrintAssignDetail(const vector< pair<int, int> >& assignment, const vector< vector<double> >& costMat)
{
	for(int i = 0; i < assignment.size(); i++)
	{
		printf("(%d, %d): %lf\n", assignment[i].first, assignment[i].second, costMat[assignment[i].first][assignment[i].second]);
	}
	
}

vector< vector<double> >  ReadWLFromFile(const char* filename, int TRX_NUM)
{
	vector< vector<double> > TRx_WLs(TRX_NUM, vector<double> (CHAN_NUM*2,0));
	FILE* fp = fopen(filename, "rt");
	if(!fp)
	{
		printf("Open file failure\n");
		return TRx_WLs;
	}
	for(int i = 0; i < TRX_NUM; i++)
	{
		for(int j = 0; j < CHAN_NUM*2 - 1; j++)
		{
			fscanf(fp, "%lf,", &TRx_WLs[i][j]);
			//printf("%lf, ", TRx_WLs[i][j]);
		}
		fscanf(fp, "%lf\n", &TRx_WLs[i][CHAN_NUM*2 - 1]);
//		printf("\n");
	}
	fclose(fp);	

	return TRx_WLs;
}

vector< vector<double> > GenerateCostMat(const vector< vector<double> >& TRx_WLs, double channel_spacing)
{
	int TRX_NUM = TRx_WLs.size();
	vector< vector<double> > costMat(TRX_NUM, vector<double> (TRX_NUM, 0));

    for(int i = 0; i < TRX_NUM; i++)
    {
        for(int j = 0; j < TRX_NUM; j++)
        {
            if(i == j) continue;
            costMat[i][j] = TuningDistTRx(TRx_WLs[i], TRx_WLs[j], channel_spacing);
        }
    }

	return costMat;
}

vector< pair<int, int>> AssignGreedy(const vector< vector<double> >& costMat)
{
	int TRX_NUM = costMat.size();
	vector< pair<int, int> > assignment (floor(TRX_NUM/2.0), pair<int, int> (-1, -1));
	pair<int, int> currPair (-1,-1);
    double min_cost = 1e8;
    vector<bool> usedFlag(TRX_NUM, false);
    for(int k = 0; k < floor(TRX_NUM/2); k++)
    {
        min_cost = 1e8;
        for(int i = 0; i < TRX_NUM; i++)
        {
            if(usedFlag[i]) continue;
            for(int j = 0; j < TRX_NUM; j++)
            {
                if(i == j || usedFlag[j]) continue;
                if(costMat[i][j] < min_cost)
                {
                    min_cost = costMat[i][j];
                    currPair.first = i;
                    currPair.second = j;
                }
            }
        }
		usedFlag[currPair.first] = true;
		usedFlag[currPair.second] = true;
		assignment[k] = currPair;
           
    //    printf("(%d, %d): %f\n", currPair.first + 1, currPair.second + 1, min_cost, channel_spacing);
    }
	return assignment;
}


int main(int argc, char** argv)
{
    //Configuration parameters
    int TRX_NUM = 31;
    double channel_spacing = 0.64*2; //0.64 nm for 80 GHz channel spacing
	double lambda1 = 0.0;
//  double wl0 = 1310; //center wavelength
//  double wl_std = 0.2;
	
	//Read the WLs
	vector< vector<double> > TRx_WLs = ReadWLFromFile("../160G_WLs.csv", TRX_NUM);

	//Generate the cost matrix
    vector< vector<double> > costMat = GenerateCostMat(TRx_WLs, channel_spacing);

	//The greedy algorithm
	//Use its result as the initial assignment for the SA algorithm
	vector< pair<int, int>> assignment = AssignGreedy(costMat);
	PrintAssignDetail(assignment, costMat);
	double initCost = AssignCost(assignment, costMat, lambda1);

	//The simulated annealing algorithm
	int pair_num = floor(TRX_NUM/2.0);
	double cost0 = initCost; //the previous-state cost
	double cost1 = 0; //the current-state cost
	
	double kT = initCost*0.5;		//Initial temperature for the SA algorithm
	double kT_scale_rate = 0.9;		//Cooling rate
	int loops_per_temp = 10000;		//Loop number for each temperature to reach equilibrium
	double E_min_per_temp0 = cost0;
	double E_min_per_temp1 = cost0;
	double E_min_global = cost0;
	vector< pair<int, int> > assignment_optimal(assignment);
	double E_change_criteria = 0.01;
	
    srand((int)time(0));
	int perturb_cnt = 0;
	while(true)
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
		if(abs(E_min_per_temp1 - E_min_per_temp0)/E_min_per_temp0 > E_change_criteria)
		{
			kT = kT*kT_scale_rate;
		}
		else break;
		
		E_min_per_temp0 = E_min_per_temp1;

	}

	printf("Initial cost: %lf\n", initCost);
	printf("Global min cost: %lf\n", E_min_global);
	printf("Optimal assignment:\n");
	PrintAssignDetail(assignment_optimal, costMat);	


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


