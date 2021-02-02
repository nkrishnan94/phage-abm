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
const int N_demes = 200; // number of demes in comiving frame
//const int N_spec = 2; // number of 'species' including phage and bacteria
const int K_bac=500; // deme size for bacteria
const int K_vir = 100; // deme size for phage - >beta*K_bac*2
float tao = 5; // lysis time in simulation steps
int beta = 100; //number of phage released with lysis
float M = .25; // Migration rate
int prof_hist = 0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 1*pow(10,5); // Run time in generations
int samp_id=0;
float alpha = 0.01;




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
    clock_t c_init = clock(); // Initial time; used to output run time
	struct tm * timeinfo;
	char buffer [80];
    time (&time_start);
	timeinfo = localtime (&time_start);
	strftime(buffer,80,"%F-%H-%M-%S", timeinfo);

	//-------------intialize data files----------------------------
    ofstream fpop, fhet, fprofp,fprofb, flog,fabs; //filenames
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
    string absName = "abs_Nb" + strKb.str()  + "_migr" + strM.str() + "_B"  +"+_tau"+strT.str()+"_alpha"+strA.str()+"_ID"+strI.str()+  termination;
    string folder = "";
    flog.open(folder+logName);
    //fhet.open(hetName);
    fpop.open(folder+velName);
    fprofp.open(folder+profpName);
    fprofb.open(folder+profbName);
    fabs.open(folder+absName);


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
    long V_deme_aux[N_demes][2] = {{0}};
    double shiftDemes = 0; // Number of demes shifted
    int shiftpop=0;
    int record_time=5000;
    vector <double> pop_hist;
    vector <double> het_hist;
    vector <double>abs_hist;
    int total_phage = int(N_demes/2)*100;
    int total_bac_no_inf=0;
    int total_bac_inf= K_bac*int(N_demes/2);
    //int tao = tao_*K_bac*beta;
    long tao_count = pow(10,int(log10(tao) + 2));
    uniform_real_distribution<double> distribution_d(0.0, 1.0);



    ///---setup iinitial population
    for(int m= 0; m<int(N_demes/2);m++){

		V_deme[m][0]=1000;



		V_deme[m][1]=1000;

    }


    ///-------Helper Funcs-----------


    //main loop
    int aux_p0;
    int aux_p1;
    int hold0;
    int hold1;
    float absAllTot;
    int acnt;

    for (unsigned int t = 0; t < N_gen; t++){
        //int acnt=0;
        //float absAllTot=0;
        
        //cout<<"hi"<<endl;
        int m=0;
        

        //first deme
        aux_p0 =V_deme[m][0]; 
        aux_p1 =V_deme[m][1]; 
        //cout<<V_deme[0][1]<<endl;
        float M0 =((1-M/2)* V_deme[m][0] + (M/2)*V_deme[m+1][0]);
        float M1 = ((1-M/2)* V_deme[m][1] + (M/2)*V_deme[m+1][1]);
        binomial_distribution<int> distribution_M0( 1,float(M0 - int(M0) ));
        binomial_distribution<int> distribution_M1( 1,float(M1 - int(M1) ));
        V_deme[m][0] = int(M0) + distribution_M0(e);
        V_deme[m][1] = int(M1) + distribution_M1(e);

        //absorption

        //cout <<alpha*(V_deme[m][0]/ (V_deme[m][0]+V_deme[m][1]))<<endl;
        //cout<<V_deme[m][0]<<endl;

        int Bempty =0;
        for (int n=0;n<K_bac;n++){
            if (B_deme[m][n]==0){
                if(B_deme[m][n]==0){
                    Bempty+=1;
                }
 
            }

        }


        if ((V_deme[m][0]+V_deme[m][1]>0)&& (Bempty>0)){
            acnt+=1;

            binomial_distribution<int> distribution_0( (V_deme[m][0] +V_deme[m][1])*Bempty ,alpha/ K_bac);
            cout<<alpha/K_bac<<endl;

            

            int atot = distribution_0(e);

            atot=min(atot,Bempty);

            absAllTot+=float(atot)/float(V_deme[m][0] +V_deme[m][1]);



            binomial_distribution<int> distribution_1( atot,float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) );
            int a0 = distribution_1(e);
            int a1 = atot-a0;
           //cout<<a0<<endl;
            int Babs = 0;
            for (int n=(K_bac - Bempty-1);n<(K_bac - Bempty+a0);n++){
                if (B_deme[m][n] ==0){
                    B_deme[m][n] = (1) * tao_count;
                    Babs+=1;

                }
            }

            for (int n=(K_bac - Bempty+a0-1);n<(K_bac -Bempty+a0+a1);n++){
                if (B_deme[m][n] ==0){
                    B_deme[m][n] = (2) * tao_count;
                    Babs+=1;

                } 
            }
            if (Babs<(a0+a1)){


                cout<<"Timestep: "<< t << " Deme: "<< m<<"Bempty: "<< Bempty<< " Bac. Index: " <<(a0+a1)<<endl;
                cout <<"NO UNINFECTED BACTERIA AVAILABLE FOR ABSORPTION" <<t <<endl;
                exit(EXIT_FAILURE);
                
            }
        }
        


        //lysis
       
        for(int nb=0; nb< K_bac;nb++){
            if (B_deme[m][nb]>0){
                B_deme[m][nb]+=1;
            }

            if ((B_deme[m][nb]% tao_count)==tao){

                int cnt=0;

                int burst_phage=(B_deme[m][nb]-tao)/tao_count-1;
                B_deme[m][nb]=-1;

                V_deme[m][burst_phage]+=beta;


            }

        }
        
        //subsequent demes
    
        for(int m=1;m<N_demes-2;m++){
            
            
            hold0 = aux_p0;
            hold1 = aux_p1;
            aux_p0 =V_deme[m][0]; 
            aux_p1 =V_deme[m][1]; 
            float M0 = ((1-M)* V_deme[m][0] + (M/2)*V_deme[m+1][0] + (M/2)*hold0);
            float M1 = ((1-M)* V_deme[m][1] + (M/2)*V_deme[m+1][1] + (M/2)*hold1);
            binomial_distribution<int> distribution_M0( 1,float(M0 - int(M0) ));
            binomial_distribution<int> distribution_M1( 1,float(M1 - int(M1) ));

            V_deme[m][0] = int(M0) + distribution_M0(e);
            V_deme[m][1] = int(M1) + distribution_M1(e);
            //cout <<alpha*(V_deme[m][0]/ (V_deme[m][0]+V_deme[m][1]))<<endl;
        //cout<<V_deme[m][0]<<endl;

            int Bempty =0;
            for (int n=0;n<K_bac;n++){
                if (B_deme[m][n]==0){
                    if(B_deme[m][n]==0){
                        Bempty+=1;
                    }
     
                }

            }

            if ((V_deme[m][0]+V_deme[m][1]>0)&& (Bempty>0)){
                binomial_distribution<int> distribution_0( (V_deme[m][0] +V_deme[m][1])*Bempty , alpha/K_bac);
                cout<<alpha/K_bac<<endl;

                acnt+=1;
                

                int atot = distribution_0(e);
                atot=min(atot,Bempty);

                absAllTot+=float(atot)/float(V_deme[m][0] +V_deme[m][1]);
                //cout<< absAllTot<<endl;
                //cout << (V_deme[m][0] +V_deme[m][1])<<" "<<Bempty<< " "<<(V_deme[m][0] +V_deme[m][1])*Bempty *(alpha/K_bac)<<endl;
                //cout<<Bempty<<" " <<atot<<endl;
                
                //cout<< float(atot)/float(V_deme[m][0] +V_deme[m][1])<<endl;

                binomial_distribution<int> distribution_1( atot,float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) );
                int a0 = distribution_1(e);
                int a1 = atot-a0;
               //cout<<a0<<endl;
                int Babs = 0;
                for (int n=(K_bac - Bempty-1);n<(K_bac - Bempty+a0);n++){
                    if (B_deme[m][n] ==0){
                        B_deme[m][n] = (1) * tao_count;
                        Babs+=1;

                    }
                }

                for (int n=(K_bac - Bempty+a0-1);n<(K_bac -Bempty+a0+a1);n++){
                    if (B_deme[m][n] ==0){
                        B_deme[m][n] = (2) * tao_count;
                        Babs+=1;

                    } 
                }
                if (Babs<(a0+a1)){


                    cout<<"Timestep: "<< t << " Deme: "<< m<<"Bempty: "<< Bempty<< " Bac. Index: " <<(a0+a1)<<endl;
                    cout <<"NO UNINFECTED BACTERIA AVAILABLE FOR ABSORPTION" <<t <<endl;
                    exit(EXIT_FAILURE);
                    
                }
            }
            


            //lysis
           
            for(int nb=0; nb< K_bac;nb++){
                if (B_deme[m][nb]>0){
                    B_deme[m][nb]+=1;
                }
                if ((B_deme[m][nb]% tao_count)==tao){

                    int cnt=0;

                    int burst_phage=(B_deme[m][nb]-tao)/tao_count-1;
                    B_deme[m][nb]=-1;

                    V_deme[m][burst_phage]+=beta;


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
                
                cout <<V_deme[m][0]<<" " << V_deme[m][1]<< " "<<m<<endl;
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



            cout<<endl;


            for(int m=0;m<N_demes;m++){
                int total_deme=0;

                for(int n=0;n<2;n++){

                    
                    total_deme+=V_deme[m][n];

                }
                cout<<total_deme <<" ";
                


            
            }

            cout<<"timestep: "<< t<<" Het: "<< calcHet(V_deme)<< "total pop "<< total_phage<<endl;
            cout<<"demes until shift"<< shift <<endl;
			pop_hist.push_back(total_phage);
            abs_hist.push_back(absAllTot/acnt);
			het_hist.push_back(calcHet(V_deme));

		}
        absAllTot=0;
        acnt=0;

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
        fabs <<dt*record_time<< " " << abs_hist[dt] <<endl;
    	fhet <<dt*record_time<< " " << het_hist[dt] <<endl;

    }


    ///final out put
    
    time_t time_end;
    clock_t c_fin = clock(); // Stop clock
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< double(c_fin - c_init)/CLOCKS_PER_SEC <<endl;

    cout << "Finished in " << double(c_fin - c_init)/CLOCKS_PER_SEC << " seconds \n";

    


    return 0;
}



