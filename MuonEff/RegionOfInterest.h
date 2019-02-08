#ifndef RegionOfInterest_h
#define RegionOfInterest_h

double ROI_func(int region, double eget){
	double p[5];

	if(region == 1){
		// fir for median
		p[0] = 0.000166057;
		p[1] = -0.396034;
		p[2] = -1.01227;
		p[3] = -0.0444885;
		p[4] = 1.34528;
	}
	if(region == 2){
		// fir for median
		p[0] = 0.000103963;
		p[1] = -0.042149;
		p[2] = -1.04393;
		p[3] = -0.097035;
		p[4] = 3.13631;
	}


	return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
}

#endif 
