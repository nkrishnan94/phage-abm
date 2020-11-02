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



//define key parameters
int tao = 15; // lysis time in simulation steps
int beta = 40; //number of phage released with lysis
const int N_demes = 300; // number of demes in comiving frame
//const int N_spec = 2; // number of 'species' including phage and bacteria
const int K_bac = 100; // deme size for bacteria
const int K_vir = 3000; // deme size for phage - >beta*K_bac
float M = .15; // Migration rate
int prof_hist = 0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 9000; // Run time in generations


//////-------HELPER FUNCTIONS ------


float calcHet(long arr [N_demes][K_vir]){

	int cnt =  0 ;
	long double H = 0.0;

	for(int m =0; m<N_demes;m++){
		float deme_pop_a1=0;
		float deme_pop_a2=0;
		for(int n=0; n<K_vir;n++){
			if (arr[m][n]==1){
				deme_pop_a1+=1;
			}
			if (arr[m][n]==1){
				deme_pop_a2+=1;
			}

		}
		if((deme_pop_a1+deme_pop_a2)>0){
			cnt+=1;
			H+=float(2*deme_pop_a1*deme_pop_a2)/ float(pow(deme_pop_a1+deme_pop_a2,2));
		}


	}

	return float(H)/float(cnt);


}

int main (int argc, char * argv[]){
	using namespace std;
	///---------read in simulation parameters from out put flags
	int c;
    while ((c = getopt (argc, argv, "B:V:M:b:t:H")) != -1)
    {

        if (c == 'b')
            beta = atof(optarg); // migration probability
        else if (c == 't')
            tao = atof(optarg); // migration probability
        else if (c == 'H')
            prof_hist = atoi(optarg); //keep track of profile through time 
    }

    //----------initialize time-------------
	time_t time_start;
	struct tm * timeinfo;
	char buffer [80];
    time (&time_start);
	timeinfo = localtime (&time_start);
	strftime(buffer,80,"%F-%H-%M-%S", timeinfo);

	//-------------intialize data files----------------------------
    ofstream fpop, fhet, fprofp,fprofb, flog; //filenames
    ostringstream strTime, strKb, strKv, strM,strDeme; //paramater strings
    strKb << K_bac;
    strKv << K_vir;
    strTime << buffer;
    strM << M;
    strDeme << N_demes;

    
    string termination = "_demes" + strDeme.str() +"_" +strTime.str() +".txt";
    string velName = "velocity_Nb" + strKb.str() + "_Nv" + strKv.str() + "_migr" + strM.str() + "_B" + termination;
    string hetName = "hetero_Nb" + strKb.str() + "_Nv" + strKv.str() + "_migr" + strM.str() + "_B"  + termination;
    string profpName = "profile_phage_Nb" + strKb.str() + "_Nv" + strKv.str() + "_migr" + strM.str()  + termination;
    string profbName = "profile_bac_Nb" + strKb.str() + "_Nv" + strKv.str() + "_migr" + strM.str()  + termination;
    string logName = "log_Nb" + strKb.str() + "_Nv" + strKv.str() + "_migr" + strM.str() + "_B"  + termination;
    flog.open(logName);
    fhet.open(hetName);
    fpop.open(velName);
    fprofp.open(profpName);
    fprofb.open(profbName);


    /// ---------------setup RNG-----------------
	//const gsl_rng_type * T;
	//gsl_rng * r;
	//gsl_rng_env_setup();
	//T = gsl_rng_mt19937;
	//r = gsl_rng_alloc(T);
	//int sysRandom;
	//gsl_rng_set(r, sysRandom);
	int seed =3133;
	int tRand = time(NULL);
	mt19937 e(tRand);

	///------additional parameters
	long B_deme[N_demes][K_bac] = {{0}};// Keep track of Bacteria -> Healthy, infected, lysed
	long V_deme[N_demes][K_vir] = {{0}};// //Keep track of N_spec species of phage
    double shiftDemes = 0; // Number of demes shifted
    int record_time=1000;
    vector <double> pop_hist;
    vector <double> het_hist;
    ///---setup iinitial population
    for(int m= 0; m<int(N_demes/2);m++){
    	for(int n =0; n< int(K_vir/4);n++){
    		V_deme[m][n]=1;

		}	
    	for(int n =int(K_vir/4); n< int(K_vir/2);n++){
    		V_deme[m][n]=2;

		}	
    }


    ///-------Helper Funcs-----------


    //main loop
    for (unsigned int t = 0; t < N_gen; t++){
        

    	//migration
	    uniform_int_distribution<int> distribution_Dm(0, N_demes);
		int r_deme = distribution_Dm(e);
		uniform_int_distribution<int> distribution_Im(0, K_vir);
		int r_ind = distribution_Im(e);
		int r_phage = V_deme[r_deme][r_ind];
		int mig_deme;

		if (r_phage>0){
			uniform_real_distribution<double> distribution_d(0.0, 1.0);
			double p_mig = distribution_d(e);
			//cout<<r_deme<<endl;


			if(r_deme>0){
				if(p_mig<M){
					double r_dir = distribution_d(e);
					int mig_deme;
					if (r_dir<.5){
						mig_deme = r_dir -1;
					} else{
						mig_deme = r_dir +1;

					}

					int found_empty =0;
					int cnt=0;
					while((found_empty==0)&&(cnt<K_vir)){
						if(V_deme[mig_deme][cnt]==0){
							found_empty=1;
							V_deme[mig_deme][cnt] = r_phage;
							V_deme[r_deme][r_ind]=0;

						} else{
							cnt+=1;
						}

					}
				}
			}

			if(r_deme==0){
				if(p_mig<M/2){

					int found_empty =0;
					int cnt=0;
					mig_deme=1;
					while((found_empty==0)&&(cnt<K_vir)){
						if(V_deme[mig_deme][cnt]==0){	
							found_empty=1;
							V_deme[mig_deme][cnt] = r_phage;
							V_deme[r_deme][r_ind]=0;

						} else{
							cnt+=1;
						}
					}

				}

			}
		}

    	//infection
	    uniform_int_distribution<int> distribution_Di(0, N_demes);
		int r_demei = distribution_Di(e);
		uniform_int_distribution<int> distribution_Ii(0, K_vir);
		int r_indi = distribution_Ii(e);
		int r_phagei = V_deme[r_demei][r_indi];
		if(r_phagei>0){
			int found_bac =0;
			int cnt=0;
			while((found_bac==0)&&(cnt<K_bac)){
				if(B_deme[r_demei][cnt]==0){
					found_bac =1;
					B_deme[r_demei][cnt]=100*r_phagei;
					V_deme[r_demei][r_ind] =0;
				} else{
					cnt+=1;
				}

			}

		}


		//lysis
		for(int m=0; m< N_demes;m++){
			for(int n=0; n< K_bac;n++){
				if (B_deme[m][n]>0){
					B_deme[m][n]+=1;
				}
				if ((B_deme[m][n]% 100)==tao){
					int found_empty=0;
					int cnt;
					int burst_phage=int((B_deme[m][n]-tao)/100);
					while(found_empty==0){
						if(V_deme[m][cnt]==0){
							found_empty=1;
							B_deme[m][n]=-1;

						} else{
							cnt+=1;
						}

					}
					int b_cnt=0;
					while((b_cnt<beta)){

                        V_deme[m][cnt] = burst_phage;
                        cnt+=1;
                        b_cnt+=1;
                        

					}

				}

			}

		}


		//shift populatoion
		int tot_pop=0;
		for(int m=0;m<N_demes;m++){
			for(int n=0;n<K_vir;n++){
				if (V_deme[m][n]>0){
					tot_pop+=1;
				}

			}

		
		}
		//cout<<tot_pop<<endl;
		//cout <<tot_pop<<endl;

		///segmentation fault here!!!
		int shift = int(tot_pop/K_vir - N_demes/2);

		if (shift>0){
            //cout<<shift<<endl;
            
			for(int s =shift; s<N_demes;s++){
				for(int p=0;p<K_vir;p++){
					V_deme[s-shift][p]=V_deme[s][p];
				}
				for(int b=0;b<K_bac;b++){
					B_deme[s-shift][b]=B_deme[s][b];
				}

			}
			for(int s = N_demes-shift; s<N_demes;s++){
				for(int p=0;p<K_vir;p++){
					V_deme[s][p]=0;
				}
				for(int b=0;b<K_bac;b++){
					B_deme[s][b]=0;
				}

			}
			shiftDemes+=shift;
		}

		if(t%record_time ==0){
            cout<< t<<" "<< shift<<endl;
			pop_hist.push_back(shiftDemes+tot_pop/K_vir);
			het_hist.push_back(calcHet(V_deme));
		}

	}
	//////-----write data files-------------
    for (int m = 0; m < N_demes; m++){
    
        int p_vac=0;
        int p_a1=0;
        int p_a2=0;
        for(int n = 0; n < K_vir; n++){
        	if (V_deme[m][n]==0){
        		p_vac+=1;
        	}
        	if (V_deme[m][n]==1){
        		p_a1+=1;
        	}
        	if (V_deme[m][n]==2){
        		p_a2+=1;
        	}	
    	}
        fprofp << m <<" " <<p_vac<<" " <<p_a1 <<" "<< p_a2 <<endl;

        int B_health=0;
        int B_inf=0;
        int B_lys=0;

        for(int n = 0; n < K_bac; n++){
        	if (B_deme[m][n]==0){
        		B_health+=1;
        	}
        	if (B_deme[m][n]>0){
        		B_inf+=1;
        	}
        	if (B_deme[m][n]==-1){
        		B_lys+=1;
        	}	
    	}
    	fprofb << B_health <<" " <<B_inf<<" " <<B_lys <<endl;



    }

    for(int t=0; t <N_gen/record_time;t++){
    	fpop <<t*record_time<< " " << pop_hist[t] <<endl;
    	fhet <<t*record_time<< " " << het_hist[t] <<endl;

    }


    ///final out put
    cout << "Finished!" << "\n";
    time_t time_end;
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< run_time <<endl;

    cout << "Finished in " << run_time << " seconds \n";




    return 0;
}




