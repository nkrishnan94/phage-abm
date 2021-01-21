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
const int N_demes = 300; // number of demes in comiving frame
//const int N_spec = 2; // number of 'species' including phage and bacteria
const int K_bac=50; // deme size for bacteria
const int K_vir = 100; // deme size for phage - >beta*K_bac*2
float tao = 100000; // lysis time in simulation steps
int beta = 50; //number of phage released with lysis
float M = .25; // Migration rate
int prof_hist = 0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 1*pow(10,7); // Run time in generations
int samp_id=0;
float alpha = 0.001;



//double B_frac = .4;


//////-------HELPER FUNCTIONS ------


float calcHet(long arr [N_demes][2]){

	double cnt =  0 ;
	long double H = 0.0;
	int deme_pop_a1;
	int deme_pop_a2;

	for(int m =0; m<N_demes;m++){
		deme_pop_a1=arr[m][0];
		deme_pop_a2=arr[m][1];

		if((deme_pop_a1+deme_pop_a2)>0){
			cnt+=1;
			H+=float(2*deme_pop_a1*deme_pop_a2)/ float(pow(deme_pop_a1+deme_pop_a2,2));
		}


	}

	return float(H)/float(cnt);


}

int main (int argc, char * argv[]){
	using namespace std;
	//coutcout <<tao_count<<endl;
	///---------read in simulation parameters from out put flags
	int c;
    while ((c = getopt (argc, argv, "b:t:i:H")) != -1)
    {

        if (c == 'b')
            beta = atof(optarg); // migration probability
        else if (c == 't')
            tao = atof(optarg); // migration probability
        else if (c == 'i')
            samp_id = atoi(optarg); //keep track of profile through time
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
    ostringstream strTime, strKb, strM,strB,strT, strDeme,strA,strI; //paramater strings
    strKb << K_bac;
    strTime << buffer;
    strB<< beta;
    strT << tao;
    strM << M;
    strI << samp_id;
    strA << alpha;
    strDeme << N_demes;

    
    string termination = "_demes" + strDeme.str() +"_" +strTime.str() +".txt";
    string velName = "velocity_Nb" + strKb.str()  + "_migr" + strM.str() +"+_tau"+strT.str()+"_alpha"+strA.str() +"_ID"+strI.str()+  termination;
    string hetName = "hetero_Nb" + strKb.str()  + "_migr" + strM.str() +"+_tau"+strT.str()+"_alpha"+strA.str() +"_ID"+strI.str()+  termination;
    string profpName = "profile_phage_Nb" + strKb.str()  + "_migr" + strM.str()  +"+_tau"+strT.str()+"_alpha"+strA.str()+"_ID"+strI.str()+ termination;
    string profbName = "profile_bac_Nb" + strKb.str()  + "_migr" + strM.str() +"+_tau"+strT.str()+"_alpha"+strA.str() +"_ID"+strI.str()+ termination;
    string logName = "log_Nb" + strKb.str()  + "_migr" + strM.str() + "_B"  +"+_tau"+strT.str()+"_alpha"+strA.str()+"_ID"+strI.str()+  termination;
    string folder = "";
    flog.open(folder+logName);
    //fhet.open(hetName);
    fpop.open(folder+velName);
    fprofp.open(folder+profpName);
    fprofb.open(folder+profbName);


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
	long V_deme[N_demes][2] = {{0}};// //Keep track of N_spec species of phage
    double shiftDemes = 0; // Number of demes shifted
    int shiftpop=0;
    int record_time=100000;
    vector <double> pop_hist;
    vector <double> het_hist;
    int total_phage = int(N_demes/2)*100;
    int total_bac_no_inf=0;
    int total_bac_inf= K_bac*int(N_demes/2);
    //int tao = tao_*K_bac*beta;
    long tao_count = pow(10,int(log10(tao) + 2));
    uniform_real_distribution<double> distribution_d(0.0, 1.0);



    ///---setup iinitial population
    for(int m= 0; m<int(N_demes/2);m++){

		V_deme[m][0]=100;



		V_deme[m][1]=100;

    }


    ///-------Helper Funcs-----------


    //main loop
    for (unsigned int t = 0; t < N_gen; t++){
        if (total_phage>1){

            //migration
            uniform_int_distribution<int> distribution_ind(0, total_phage-1);
            int phage_ind = distribution_ind(e)+1;
            int phage_cnt=0;
            int phage_found=0;
            int r_deme;
            int r_phage;
            for(int m =0;m<N_demes;m++){
                if (phage_found==0){

                    phage_cnt+=V_deme[m][0];
                    if (phage_cnt+1>phage_ind){
                        r_deme=m;
                        r_phage = 0;
                        phage_found=1;

                    }
                }

                if (phage_found==0){
                    phage_cnt+=V_deme[m][1];

                    if (phage_cnt+1>phage_ind){
                        r_deme=m;
                        r_phage = 1;
                        phage_found=1;

                    }

                }
                

            }

            //cout<< r_phage<<endl;
            //cout<<total_phage<<endl;
            //cout <<phage_ind<<endl;
            //cout <<phage_cnt<<endl;
            //cout<<phage_found<<endl;
            //cout <<r_deme<<endl;
            //int r_phage = V_deme[r_deme][r_ind];
            int mig_deme;
            
            //uniform_real_distribution<double> distribution_d(0.0, 1.0);
            double p_mig = distribution_d(e);
            //cout<<r_deme<<endl;


            if(r_deme>0){
                if(p_mig<M){
                    double r_dir = distribution_d(e);
                    if (r_dir<.5){
                        mig_deme = r_deme - 1;
                    } else{
                        mig_deme = r_deme + 1;

                    }
                    if(V_deme[r_deme][r_phage]<1){
                        cout<<"migration error"<<endl;
                        
                    }
                    V_deme[r_deme][r_phage]-=1;
                    V_deme[mig_deme][r_phage]+=1;


                }
            }

            if(r_deme==0){
                //cout<<"r_deme =0 "<<endl;
                if(V_deme[r_deme][r_phage]<1){
                    cout<<"migration error"<<endl;
                    cout<< V_deme[r_deme][0]<<" "<<V_deme[r_deme][1]<<endl;
                    cout<< r_deme<< " "<<r_phage<<endl;
                    cout<<phage_ind<<" "<<phage_cnt<<endl;
                    
                }
                if(p_mig<M/2){

                    mig_deme=1;

                    V_deme[0][r_phage]-=1;
                    V_deme[1][r_phage]+=1;


                }



            }
            

            //infection
            //uniform_int_distribution<int> distribution_indi(0, total_phage-1);
            int phage_indi = distribution_ind(e) +1;
            //uniform_real_distribution<double> distribution_d(0.0, 1.0);
            double p_inf = distribution_d(e);
            int phage_cnti=0;
            int phage_foundi=0;
            int r_demei;
            int r_phagei;
            if (p_inf<alpha){
                for(int m =0;m<N_demes;m++){
                    if (phage_foundi==0){

                        phage_cnti+=V_deme[m][0];
                        if (phage_cnti+1>phage_indi){
                            r_demei=m;
                            r_phagei = 0;
                            phage_foundi=1;

                        }
                    }
                    if (phage_foundi==0){
                        phage_cnti+=V_deme[m][1];

                        if (phage_cnti+1>phage_indi){
                            r_demei=m;
                            r_phagei = 1;
                            phage_foundi=1;

                        }

                    }
                }
                //cout<<r_phagei<<endl;


                int B_found=0;
                int bac_ind;
                for (int n=0;n<K_bac;n++){
                    if (B_found==0){
                        if(B_deme[r_demei][n]==0){

                            B_found=1;
                            bac_ind=n;
                        }


                    }

                }

                if (B_found==1){
                    B_deme[r_demei][bac_ind] = (1+r_phagei)*tao_count;
                    if(V_deme[r_demei][r_phagei]<1){
                        cout<<"infection error"<<endl;
                        cout<< V_deme[r_demei][0]<<" "<<V_deme[r_demei][1]<<endl;
                        cout<< r_demei<< " "<<r_phagei<<endl;
                        cout<<phage_indi<<" "<< total_phage<<" "<<phage_cnti<<endl;
                        
                        
                    }
                        
                    V_deme[r_demei][r_phagei]-=1;

                }
            }
        }







            //lysis
        for(int m=0; m< N_demes;m++){
            for(int nb=0; nb< K_bac;nb++){
                if (B_deme[m][nb]>0){
                    B_deme[m][nb]+=1;
                }
                if ((B_deme[m][nb]% tao_count)==tao){
                    //if (m==10){
                     //   cout <<V_deme[m][0]<<" "<<V_deme[m][1]<<endl;
                   // }
                    //cout<<B_deme[m][nb]<<endl;
                    
                    //int found_empty=0;
                    int cnt=0;

                    int burst_phage=(B_deme[m][nb]-tao)/tao_count-1;
                    //cout << burst_phage<<endl;
                    B_deme[m][nb]=-1;
                    //cout<<m<<endl;
                    V_deme[m][burst_phage]+=beta;
                    //if (m==10){
                     //   cout <<V_deme[m][0]<<" "<<V_deme[m][1]<<" "<<burst_phage<<endl;
                    //}
                    //total_phage+=beta;

                }

            }

        }
        

		//shift populatoion
		total_phage=0;
		int last_deme=0;
		for(int m=0;m<N_demes;m++){

			total_phage+=V_deme[m][0]+V_deme[m][1];
			if (V_deme[m][0]+V_deme[m][1]>0){
				last_deme=m;


			}
            if(( V_deme[m][0]<0 )|| (V_deme[m][1]<0)){
                cout <<"deme population negative, timestep " <<t <<endl;
                
                //cout <<V_deme[m][0]<<" " << V_deme[m][1]<<endl;
                exit(EXIT_FAILURE);
            }


		
		}
		//cout<<tot_pop<<endl;
		//cout <<tot_pop<<endl;

		///segmentation fault here!!!
        int shift = 0;
        //if (last_deme >= N_demes*(3/4)){
        shift = last_deme-10-int(N_demes/2);
            
        
		//cout<<last_deme<<endl;


		if (shift>0){
            //cout<<shift<<endl;
            //shiftpop+= V_deme[s][0]+V_deme[s][1];
            for(int s =0;s<shift;s++){

            	shiftpop+= V_deme[s][0]+V_deme[s][1];

        	}
            
			for(int s =shift; s<N_demes;s++){
				
				V_deme[s-shift][0]=V_deme[s][0];
				V_deme[s-shift][1]=V_deme[s][1];
				
				for(int b=0;b<K_bac;b++){
					B_deme[s-shift][b]=B_deme[s][b];
				}

			}
			for(int s = N_demes-shift; s<N_demes;s++){
				V_deme[s][0]=0;
				V_deme[s][1]=0;

				for(int b=0;b<K_bac;b++){
					B_deme[s][b]=0;
				}

			}
			
		}
        //total_phage-=shiftpop;



		if(t%record_time ==0){
            for(int m=0;m<N_demes;m++){
                int total_deme=0;

                for(int n=0;n<2;n++){

                    
                    total_deme+=V_deme[m][n];

                }
                cout<<total_deme <<" ";

            
            }
            cout<<"timestep: "<< t<<" Het: "<< calcHet(V_deme)<< "total pop: "<< total_phage<< "demes until shift: "<< shift<<endl;
			pop_hist.push_back(shiftpop+total_phage);
			//het_hist.push_back(calcHet(V_deme));
		}

	}
	
	//////-----write data files-------------
    for (int m = 0; m < N_demes; m++){
    

        fprofp << m <<" "<<V_deme[m][0]<<" " <<V_deme[m][1] <<endl;

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
    	fprofb <<m<<" "<< B_health <<" " <<B_inf<<" " <<B_lys <<endl;



    }
    //std::cout <<"hi"<<std::endl;
    for(int dt=0; dt <int(N_gen/record_time);dt++){

    	fpop <<dt*record_time<< " " << pop_hist[dt] <<endl;
    	//fhet <<t*record_time<< " " << het_hist[t] <<endl;

    }


    ///final out put
    
    time_t time_end;
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< run_time <<endl;

    cout << "Finished in " << run_time << " seconds \n";

    


    return 0;
}



