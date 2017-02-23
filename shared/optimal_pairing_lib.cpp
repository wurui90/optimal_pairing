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
#include <boost/random/normal_distribution.hpp>
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
			ret[i][j] = GV + j*channel_spacing + LV_gen(); //Tx WLs
			ret[i][j + channel_num] = GV + TRxV + j*channel_spacing + LV_gen(); //Rx_WLs
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
void PrintAssignDetail(const vector< pair<int, int> >& assignment, const vector< vector<double> >& costMat)
