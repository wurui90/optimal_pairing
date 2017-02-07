#include <stdio.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <utility>

using namespace std;

bool myLessThan (float i, float j) {return (i < j);}

double dist(vector<double> wl1, vector<double> wl2)
{
    double ret = 0;
    for(int i = 0; i < wl1.size(); i++)
    {
        ret += abs(wl1[i] - wl2[i]);
    }
    return ret;
}

pair<int, int> greedyFindPair(const vector< vector<double> >& Tx, const vector< vector<double> >& Rx, vector<bool>& TxUsedFlag, vector<bool> RxUsedFlag)
{
    int TRx_num = Tx.size();
    vector< vector<double> > costMat(TRx_num, vector<double> (TRx_num, 0));

    double min_cost = 1e8;
    pair<int, int> min_index;
    for(int i = 0; i < TRx_num; i++)
    {
        if(TxUsedFlag[i]) continue;
        for(int j = 0; j < TRx_num; j++)
        {
            if(RxUsedFlag[j]) continue;
            costMat[i][j] = dist(Tx[i], Rx[j]);
            if(costMat[i][j] < min_cost)
            {
                min_cost = costMat[i][j];
                min_index.first = i;
                min_index.second = j;
            }
        }
    }

    TxUsedFlag[min_index.first] = true;
    RxUsedFlag[min_index.second] = true;
    return min_index;
}


//main func args
//1: TRx number
//2: WDM channel number
int main(int argc, char** argv)
{
    int TRX_NUM = 100;
    int CHAN_NUM = 2;

    TRX_NUM = atoi(argv[1]);
    CHAN_NUM = atoi(argv[2]);

    //Randomly generate the WLs
    printf("Synthetic WLs:\n");
    srand((int)time(0));
    vector< vector<double> > Tx_WLs(TRX_NUM, vector<double> (CHAN_NUM,0));
    vector< vector<double> > Rx_WLs(TRX_NUM, vector<double> (CHAN_NUM,0));
    vector<bool> TxUsedFlag(TRX_NUM, false);
    vector<bool> RxUsedFlag(TRX_NUM, false);

    for(int i = 0; i < TRX_NUM; i++)
    {
        printf("Device#%d:\n", i);
        printf("Tx:");
        for(int j = 0; j < CHAN_NUM; j++)
        {
            Tx_WLs[i][j] = rand()%100;
            printf("%f\t", Tx_WLs[i][j]);
        }
        printf("\n");
        printf("Rx:");
        for(int j = 0; j < CHAN_NUM; j++)
        {
            Rx_WLs[i][j] = rand()%100;
            printf("%f\t", Rx_WLs[i][j]);
        }

        printf("\n");
    }

    //Iteratively find the TRx pair with the smallest cost
    pair<int, int> currPairIndex;
    for(int i = 0; i < TRX_NUM; i++)
    {
        currPairIndex = greedyFindPair(Tx_WLs, Rx_WLs, TxUsedFlag, RxUsedFlag);
        printf("(%d, %d): %f\n", currPairIndex.first, currPairIndex.second, dist(Tx_WLs[currPairIndex.first], Rx_WLs[currPairIndex.second]));
    }

    return 0;	
}
