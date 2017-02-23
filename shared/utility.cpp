#include <vector>
#include <numeric>

using namespace std;

double Vmean(vector<double> v)
{
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = sum / v.size();
	return mean;
}

double Vstd(vector<double> v)
{
double mean = Vmean(v);
double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
return stdev;
}
