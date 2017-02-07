#include <stdio.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#define RANGE 100

using namespace std;

bool myLessThan (float i, float j) {return (i < j);}

int main()
{
    srand((int)time(0));
	vector<float> v(10, 0);

    for(int i = 0; i < v.size(); i++)
    {
        v[i] = rand()%RANGE;
        printf("%f\n", v[i]);
    }

    sort(v.begin(), v.end(), myLessThan);
    for(int i = 0; i < v.size(); i++)
    {
        printf("%f\n", v[i]);
    }

	return 0;	
}
