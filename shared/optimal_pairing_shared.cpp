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
#include "../shared/utility.cpp"

#define CHAN_NUM 5
#define PRINT_DB(X) printf("X: %lf\n", X);

using namespace std;

//% Generate ramdon WLs based on the GV and LV variation model
// Return a vector of 2*channel_num, Tx_WLs array and Rx_WLs array
vector< vector<double>> gen_syn_WLs_model(double GV_std, double LV_std, double TRxV_std, int channel_num, double channel_spacing, double base_wl, int device_num)
{
	vector< vector<double>> ret(device_num, vector<double> (channel_num*2, 0.0));
	boost::mt19937 *rng = new boost::mt19937();
	rng->seed(time(NULL));

	boost::normal_distribution<> GV_dist(0.0, GV_std);
	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > GV_gen(*rng, GV_dist);

	boost::normal_distribution<> LV_dist(0.0, LV_std);
	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > LV_gen(*rng, LV_dist);

	boost::normal_distribution<> TRxV_dist(0.0, TRxV_std);
	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > TRxV_gen(*rng, TRxV_dist);

	for(int i = 0; i < device_num; i++)
	{
		double GV = GV_gen();
		double TRxV = TRxV_gen();
		for(int j = 0; j < channel_num; j++)
		{
			ret[i][j] = GV + j*channel_spacing + LV_gen() + base_wl; //Tx WLs
			ret[i][j + channel_num] = GV + TRxV + j*channel_spacing + LV_gen() + base_wl; //Rx_WLs
		}
	}

	return ret;
}



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
		if(assignment[i].first < 0 || assignment[i].second < 0) {pair_num--;  continue;}
		A_sum += costMat[assignment[i].first][assignment[i].second];
		A_sqr_sum += costMat[assignment[i].first][assignment[i].second]*costMat[assignment[i].first][assignment[i].second];
	}
	return A_sum/pair_num + lambda1*sqrt(A_sqr_sum/pair_num - A_sum*A_sum/pair_num/pair_num);
	
}

// Print the details for a given assignment, index starts from 1
// Return the tuning cost array
vector<double> PrintAssignDetail(const vector< pair<int, int> >& assignment, const vector< vector<double> >& costMat)
{
	vector<double> ret;
	for(int i = 0; i < assignment.size(); i++)
	{
		if(assignment[i].first < 0 || assignment[i].second < 0)
		{
			printf("(%d, %d): 0\n", assignment[i].first, assignment[i].second);
		}
		else
		{
		printf("(%d, %d): %lf\n", assignment[i].first + 1, assignment[i].second + 1, costMat[assignment[i].first][assignment[i].second]);
		ret.push_back(costMat[assignment[i].first][assignment[i].second]);
		}
	}
	printf("Mean: %f, std: %f\n", Vmean(ret), Vstd(ret));
	return ret;
}

vector< vector<double> >  ReadWLFromFile(const char* filename, int TRX_NUM)
{
	vector< vector<double> > TRx_WLs(TRX_NUM, vector<double> (CHAN_NUM*2,0));
	FILE* fp = fopen(filename, "rt");
	if(!fp)
	{
		printf("Open file failure\n");
		return TRx_WLs;
		exit(-2);
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
	vector< pair<int, int> > assignment (ceil(TRX_NUM/2.0), pair<int, int> (-1, -1));
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
	//handle the odd TRx case
	if(TRX_NUM%2 == 1)
	{
		for(int i = 0; i < TRX_NUM; i++)
		{
			if(!usedFlag[i])
			{
				assignment.back().first = i;
				break;
			}
		}
	}
	return assignment;
}

