//Greedy algorithm

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <utility>
#include <random>
#include <assert.h>

#define CHAN_NUM 5

using namespace std;

bool myLessThan (float i, float j) {return (i < j);}

//tuning cost of aligning wl1 and wl2 with given wavelength grid
double tuning_dist(vector<double> wl1, vector<double> wl2, double channel_spacing)
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
double dist(vector<double> wl1, vector<double> wl2, double channel_spacing)
{
   // return abs(wl1[0] - wl2[0]);
    return tuning_dist(vector<double> (wl1.begin(), wl1.begin() + CHAN_NUM), vector<double> (wl2.begin() + CHAN_NUM, wl2.end()), channel_spacing) + tuning_dist(vector<double> (wl1.begin() + CHAN_NUM, wl1.end()), vector<double> (wl2.begin(), wl2.begin() + CHAN_NUM), channel_spacing);
}

pair<int, int> greedyFindPair(vector< vector<double> > TRx, vector<bool> &usedFlag, double channel_spacing)
{
    int TRx_num = TRx.size();
    vector< vector<double> > costMat(TRx_num, vector<double> (TRx_num, 0));

    double min_cost = 1e8;
    pair<int, int> min_index;
    for(int i = 0; i < TRx_num; i++)
    {
        if(usedFlag[i]) continue;
        for(int j = 0; j < TRx_num; j++)
        {
            if(usedFlag[j]) continue;
            costMat[i][j] = dist(TRx[i], TRx[j], channel_spacing);
            if(i != j && costMat[i][j] < min_cost)
            {
                min_cost = costMat[i][j];
                min_index.first = i;
                min_index.second = j;
            }
        }
    }

    usedFlag[min_index.first] = true;
    usedFlag[min_index.second] = true;
    return min_index;
}

int main(int argc, char** argv)
{
    //Configuration parameters
    int TRX_NUM = 40;
    double channel_spacing = 0.64;
    double wl0 = 1310; //center wavelength
    double wl_std = 0.2;
	
	vector< vector<double> > TRx_WLs(TRX_NUM, vector<double> (CHAN_NUM*2,0));
//    TRX_NUM = atoi(argv[1]);
   // CHAN_NUM = atoi(argv[2]);

    //Randomly generate the WLs
/*    normal_distribution<double> norm_dist(0.0, wl_std);
    default_random_engine rng;

    printf("Synthetic WLs:\n");
    srand((int)time(0));

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

	//Read the WLs
	FILE* fp = fopen("../80G_WLs.csv", "rt");
	for(int i = 0; i < TRX_NUM; i++)
	{
		for(int j = 0; j < CHAN_NUM*2 - 1; j++)
		{
			fscanf(fp, "%lf,", &TRx_WLs[i][j]);
			printf("%lf, ", TRx_WLs[i][j]);
		}
		fscanf(fp, "%lf\n", &TRx_WLs[i][CHAN_NUM*2 - 1]);
		printf("\n");
	}
	fclose(fp);	


    //Iteratively find the TRx pair with the smallest cost
    vector< vector<double> > costMat(TRX_NUM, vector<double> (TRX_NUM, 0));
    vector<bool> usedFlag(TRX_NUM, false);

    //Generate the cost matrix
    for(int i = 0; i < TRX_NUM; i++)
    {
        for(int j = 0; j < TRX_NUM; j++)
        {
            if(i == j) continue;
            costMat[i][j] = dist(TRx_WLs[i], TRx_WLs[j], channel_spacing);
        }
    }

    pair<int, int> currPairIndex (-1,-1);
    double min_cost = 1e8;
    for(int k = 0; k < TRX_NUM/2; k++)
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
                    currPairIndex.first = i;
                    currPairIndex.second = j;
                }
            }
        }
    usedFlag[currPairIndex.first] = true;
    usedFlag[currPairIndex.second] = true;
           
        printf("(%d, %d): %f\n", currPairIndex.first + 1, currPairIndex.second + 1, min_cost, channel_spacing);
    }

    return 0;	
}
