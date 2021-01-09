#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime> 

unsigned int NGen = 1*pow(10,6);
const int Nmax = 1*pow(10,6);
const int Nallele= 16;
float drate = .3;
float mu = .01;





int main(){
	//------
	long fitLand[Nallele] = {0};
	int max=0;

	for(int n =0; n<Nallele; n++){
        uniform_real_distribution<double> distribution_d(0.0, 1.0);
        fitLand[n] = distribution_d(e);
        if(max<fitLand[n]){

        	max = fitLand[n];
        }


	}
	for(int n =0; n<Nallele; n++){
        fitLand[n]=fitLand[n]/max;

	}
	int transMat[Nallele][Nallele] ={{0}};
	for(int nn=0; nn<Nallele; nn++){
		for(int mm=0; mm<Nallele; mm++){

			tra



	}

	//---initialize rng----
	int tRand = time(NULL);
	mt19937 e(tRand);









	return 0;
}