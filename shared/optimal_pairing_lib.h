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
