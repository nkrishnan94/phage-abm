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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int main(){

	
	using namespace std;
	vector <int>  a(6,0);

	random_shuffle(begin(a), end(a));

	for(int i =0; i <6; i++){

		cout<<a[i];
	}


	return 0;
}