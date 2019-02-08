#ifndef withEM_SingleCrys_h
#define withEM_SingleCrys_h

double SW_func2_dphi_v2(int region, double eget){
 double p[5] = {};

 if(region == 1){
   p[0] = -0.000163799;
   p[1] = 0.400094;
   p[2] = -0.938122;
   p[3] = 0.0346451;
   p[4] = 1.33086;
 }

 if(region == 2){
   p[0] = 0.000383048;
   p[1] = 0.179742;
   p[2] = -1.11884;
   p[3] = -0.548157;
   p[4] = 1.30843;
 }


 return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);
}

double SW_func2_deta_v2(double eget){
/*
 double p[1] = {};

 p[0] = -1.67124e-05;

 return p[0]*pow(eget,0);
*/
	double p[2] = {};

	p[0] = 0.00106159;
	p[1] = 0.00288571;

	return p[0]*pow(eget,0) + p[1]*pow(eget,0)*exp(-pow(eget,0));
}


#endif

