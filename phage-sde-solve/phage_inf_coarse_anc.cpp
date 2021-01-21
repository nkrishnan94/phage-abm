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
const int N_demes = 100; // number of demes in comiving frame
//const int N_spec = 2; // number of 'species' including phage and bacteria
const int K_bac=100; // deme size for bacteria
const int K_vir = 100; // deme size for phage - >beta*K_bac*2
float tao = 200; // lysis time in simulation steps
int beta = 50; //number of phage released with lysis
float M = .25; // Migration rate
int prof_hist = 0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 1*pow(10,3); // Run time in generations
int samp_id=0;
float alpha = 0.05;
unsigned int assign_gene_time = 1*pow(10,4); // Run time in generations



//double B_frac = .4;


//////-------HELPER FUNCTIONS ------


/*float calcHet(long arr [N_demes][N_demes]){

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


}*/


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
    ofstream fpop, fhet, fprofp,fprofb, flog,fmark; //filenames
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
    string markName = "markers_" + strKb.str()  + "_migr" + strM.str() + "_B" +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str()+ "_ID"+strI.str()+termination;
    string folder = "";
    flog.open(folder+logName);
    //fhet.open(hetName);
    fpop.open(folder+velName);
    fprofp.open(folder+profpName);
    fprofb.open(folder+profbName);
    fmark.open(folder+markName);


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
	long V_deme[N_demes][N_demes] = {{0}};// //Keep track of N_spec species of phage
    long V_deme_aux[N_demes][N_demes] = {{0}};
    long phage_abs_count[N_demes] = {0};
    double shiftDemes = 0; // Number of demes shifted
    int shiftpop=0;
    int record_time=500;
    vector <double> pop_hist;
    vector <double> mark_hist;
    //vector <double> het_hist;
    int total_phage = int(N_demes/2)*100;
    int total_bac_no_inf=0;
    int total_bac_inf= K_bac*int(N_demes/2);
    //int tao = tao_*K_bac*beta;
    long tao_count = pow(10,int(log10(tao) + 2));
    int t=0;
    int demes_alive=0;
    int fixed_marker;
    int N_mark = 2;
    uniform_real_distribution<double> distribution_d(0.0, 1.0);



    ///---setup iinitial population
    for(int m= 0; m<int(N_demes/2);m++){

		V_deme[m][0]=100;



		V_deme[m][1]=100;

    }


    ///-------Helper Funcs-----------


    //main loop
    int aux_p0;
    int aux_p1;
    int hold0;
    int hold1;

    while((demes_alive!=1) ||(t<assign_gene_time) ){
    //for (unsigned int t = 0; t < N_gen; t++){
        //cout<<"hi"<<endl;
      
        int deme_tot=0;

        for(int n=0;n<N_mark;n++){
            hold0= V_deme_aux[0][n];
            V_deme_aux[0][n] = V_deme[0][n];


            V_deme[0][n] = int((1-M/2)* V_deme[0][n] + (M/2)*V_deme[1][n]);
            deme_tot+=V_deme[0][n];

            

        }


        //absorption

        //cout <<alpha*(V_deme[m][0]/ (V_deme[m][0]+V_deme[m][1]))<<endl;
        //cout<<V_deme[m][0]<<endl;
        if ((deme_tot)>0){


            int abs_count_tot=0;
            for(int n=0; n <N_mark;n++){
                binomial_distribution<int> distribution( V_deme[0][n],alpha*float(V_deme[0][n]/ float(deme_tot) ) );

                //int a= distribution(e);
                phage_abs_count[n] = distribution(e);

                abs_count_tot+=phage_abs_count[n];

                //cout<<a0<<endl;
            }
                
            int bac_ind;

            vector <int> phage_abs0(abs_count_tot,0);
            if (phage_abs0.size()>0){
                //cout<<phage_abs0.size()<<endl;
                int last_mark_count=0;
                int mark_count=0;
                for(int n=0; n <N_mark;n++){
                    mark_count += phage_abs_count[n];
                    

                    for(int v=last_mark_count; v<(mark_count);v++){
                        phage_abs0[v]=n;

                    }
                    
                    last_mark_count = mark_count;

                //cout<<a0<<endl;
                }





                random_shuffle(begin(phage_abs0), end(phage_abs0));
                int r_phagei;

                

                for(int v=0; v<abs_count_tot;v++){
                    r_phagei = phage_abs0[v];
                    int B_found=0;
                    for (int n=0;n<K_bac;n++){
                        if (B_found==0){
                            if(B_deme[0][n]==0){
                                B_found=1;
                                bac_ind=n;
                            }

                        }

                    }

                    if ((B_found==1)&&(V_deme[0][r_phagei]>0)){
                        B_deme[0][bac_ind] = (1 + r_phagei) * tao_count;
                            
                        V_deme[0][r_phagei]-=1;


                    }
                }

            }
        }
        


        //lysis
       
        for(int nb=0; nb< K_bac;nb++){
            if (B_deme[0][nb]>0){
                B_deme[0][nb]+=1;
            }

            if ((B_deme[0][nb]% tao_count)==tao){

                int cnt=0;

                int burst_phage=(B_deme[0][nb]-tao)/tao_count-1;
                B_deme[0][nb]=-1;

                V_deme[0][burst_phage]+=beta;


            }

        }
        
        //subsequent demes
    
        for(int m=1;m<N_demes-2;m++){
            
           // hold0 = aux_p0;
            //hold1 = aux_p1;
            //aux_p0 =V_deme[m][0]; 
            //aux_p1 =V_deme[m][1]; 



     
            int deme_tot=0;

            for(int n=0;n<N_mark;n++){
                //hold0 = V_deme_aux[m-1][n];
                V_deme_aux[m][n] = V_deme[m][n];


                V_deme[m][n] = int((1-M)* V_deme[m][n] + (M/2)*V_deme[m+1][n] + (M/2)*V_deme_aux[m-1][n] );

                deme_tot+=V_deme[m][n];

                

            }


            //absorption

            //cout <<alpha*(V_deme[m][0]/ (V_deme[m][0]+V_deme[m][1]))<<endl;
            //cout<<V_deme[m][0]<<endl;
            if ((deme_tot)>0){


                int abs_count_tot=0;
                for(int n=0; n <N_mark;n++){
                    binomial_distribution<int> distribution( V_deme[m][n],alpha*float(V_deme[m][n]/ float(deme_tot) ) );

                    //int a= distribution(e);
                    phage_abs_count[n] = distribution(e);
                    abs_count_tot+=phage_abs_count[n];



                    //cout<<a0<<endl;
                }
                    
                int bac_ind;

                vector <int> phage_abs(abs_count_tot,0);
                if (phage_abs.size()>0){
                    int last_mark_count=0;

                    int mark_count=0;
                    for(int n=0; n <N_mark;n++){
                        mark_count += phage_abs_count[n];
                        //cout<<mark_count<<endl;
    

                        for(int v=last_mark_count; v<(mark_count);v++){
                            phage_abs[v]=n;

                            

                        }
                    
                        last_mark_count = mark_count;

                    //cout<<a0<<endl;
                    }
                    //cout<<mark_count<<endl;






                    random_shuffle(begin(phage_abs), end(phage_abs));
                    int r_phagei;

                    

                    for(int v=0; v<abs_count_tot;v++){
                        //if (phage_abs[v]==0){
                         //   cout<< phage_abs[v] <<endl;


                        //}
 

                        r_phagei = phage_abs[v];
                        //cout<<r_phagei<<endl;
                        int B_found=0;
                        for (int n=0;n<K_bac;n++){
                            if (B_found==0){
                                if(B_deme[m][n]==0){
                                    B_found=1;
                                    bac_ind=n;
                                }

                            }

                        }

                        if ((B_found==1)&&(V_deme[m][r_phagei]>0)){
                            B_deme[m][bac_ind] = (1 + r_phagei) * tao_count;

                                
                            V_deme[m][r_phagei]-=1;



                        }
                    }

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
                    //cout<<m<<endl;
                    B_deme[m][nb]=-1;

                    V_deme[m][burst_phage]+=beta;


                }

            }



        }

        total_phage=0;

        int last_deme=0;

        //for(int n=0;n<N_demes;n++){
        //    V_tot[n]=0;
        //}

        for(int m=0;m<N_demes;m++){
            int total_deme=0;
            int max_deme=0;
            for(int n=0;n<N_mark;n++){

                
                total_deme+=V_deme[m][n];
                //V_tot[n]+=V_deme[m][n];

                if (V_deme[m][n]<0){
                   cout <<"deme population negative, timestep " <<t <<endl;
                    //cout <<V_deme[m][0]<<" " << V_deme[m][1]<<endl;
                    exit(EXIT_FAILURE);

                }

            }
            //V_prof[m] = total_deme;

            total_phage+=total_deme;
            
            if (total_deme>0){
                last_deme=m;

            }


        
        }

        if (t>assign_gene_time){
            
            demes_alive=0;
            
            for(int m=0;m<N_demes;m++){
                int org_tot =0;

                for(int n=0;n<N_mark;n++){


                    org_tot+=V_deme[n][m];

                }
                if(org_tot>0){

                    demes_alive+=1;
                }



            }


            if (demes_alive==0){
                cout <<"NO PHAGE ORIGINS RECORDED"<< endl;

            }
            if (demes_alive==1){

                int fixed_found=0;
                for(int m=0;m<N_demes;m++){
                    int org_tot =0;

                    for(int n=0;n<N_mark;n++){

                        org_tot+=V_deme[n][m];

                    }
                    if(org_tot>0){
                        if (fixed_found==1){
                            cout<<"ONLY ONE GENE ORIGIN SHOULD BE FOUND"<<endl;
                        }
                        fixed_marker=m;
                        cout << "Phage from deme: " << m <<" FIXED" << endl;
                        fixed_found=1;
                        

                    }

                }


            }



        }
        

        
        

		//shift populatoion
		int shift = 0;
        //if (last_deme >= N_demes*(3/4)){
        shift = last_deme-10-int(N_demes/2);
            
        
        //cout<<last_deme<<endl;
        

        if (shift>0){


            //cout<<"Shift: " <<shift<<endl;
            //shiftpop+= V_deme[s][0]+V_deme[s][1];
            for(int s =0;s<shift;s++){
                for(int m=0;m<N_demes;m++){

                    shiftpop+= V_deme[s][m];



                }


            }
            
            for(int s =shift; s<N_demes;s++){
                 for(int m=0;m<N_demes;m++){


                
                    V_deme[s-shift][m]=V_deme[s][m];

                }

                
                for(int b=0;b<K_bac;b++){
                    B_deme[s-shift][b]=B_deme[s][b];
                }

            }
            for(int s = N_demes-shift; s<N_demes;s++){
                for(int m=0;m<N_demes;m++){
                    V_deme[s][m]=0;

                }
                //V_deme[s][0]=0;
                //V_deme[s][1]=0;

                for(int b=0;b<K_bac;b++){
                    B_deme[s][b]=0;
                }

            }
            
        }
        //total_phage-=shiftpop;



		if(t%record_time ==0){
            //for(int n=0;n<N_mark;n++){
            for(int m=0;m<N_demes;m++){
                int deme_tot =0;
                //int org_tot=0;

                for(int n=0;n<N_mark;n++){
                //for(int m=0;m<N_demes;m++){


                    //org_tot+=V_deme[n][m];
                    deme_tot+=V_deme[m][n];

                }
                //cout<< org_tot <<" ";
                cout<< deme_tot <<" ";
                        //for(int n=0;n<N_demes;n++){
        //    V_tot[n]=0;
        //}


            }
            cout<<endl;

            cout<<"timestep: "<< t<< "total pop "<< total_phage<<endl;
            cout<<"demes until shift"<< shift <<endl;

			pop_hist.push_back(shiftpop+total_phage);
            mark_hist.push_back(demes_alive);
			//het_hist.push_back(calcHet(V_deme));
		}

        if (t>assign_gene_time){
            
            demes_alive=0;
            
            //for(int m=0;m<N_demes;m++){
            for(int n=0;n<N_mark;n++){
                int org_tot =0;

                //for(int n=0;n<N_mark;n++){
                for(int m=0;m<N_demes;m++){


                    org_tot+=V_deme[m][n];

                }
                if(org_tot>0){

                    demes_alive+=1;
                }



            }


            if (demes_alive==0){
                cout <<"NO PHAGE ORIGINS RECORDED"<< endl;
            }

            if (demes_alive==1){

                int fixed_found=0;
                //for(int m=0;m<N_demes;m++){
                for(int n=0;n<N_mark;n++){
                    int org_tot =0;

                    for(int m=0;m<N_demes;m++){

                        org_tot+=V_deme[m][n];

                    }
                    if(org_tot>0){
                        if (fixed_found==1){
                            cout<<"ONLY ONE GENE ORIGIN SHOULD BE FOUND"<<endl;
                        }
                        fixed_marker=n;
                        cout << "Phage from deme: " << n <<" FIXED" << endl;
                        fixed_found=1;
                        

                    }

                }


            }
        }

        if(t==assign_gene_time){
            N_mark =N_demes;
            //save 'initial ' profile'

            cout<<"Saving profiles\n";
            for (int m = 0; m < N_demes; m++){
                fprofp << m <<" "<<V_deme[m][0]<<endl;


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
            //assign label to free phage
            int deme_pop=0;
            for(int m=0;m<N_demes;m++){
                int deme_pop=0;
                for(int n=0;n<N_mark;n++){
                    deme_pop+=V_deme[m][n];

                }

                for(int n=0;n<N_mark;n++){
                    V_deme[m][n]=0;

                }
                V_deme[m][m]= deme_pop;
                cout<<V_deme[m][m]<<endl;



            }
            //assign label  to absorbed phage
            for(int m=0; m< N_demes;m++){
                for(int nb=0; nb< K_bac;nb++){
                    if (B_deme[m][nb]>0){
                        B_deme[m][nb] = B_deme[m][nb]*m;


                    }

                }

            }




        }
        t+=1;

	}
	
	//////-----write data files-------------
    int deme_pop=0;
    for(int n=0;n<N_mark;n++){
    //for(int m=0;m<N_demes;m++){
        int deme_pop=0;
        //for(int n=0;n<N_mark;n++){
        for(int m=0;m<N_demes;m++){
            deme_pop+=V_deme[m][n];

        }
        cout<<deme_pop<<" ";



    }
    cout<<endl;
    
    //std::cout <<"hi"<<std::endl;
   for(int dt=0; dt <int(t/record_time);dt++){

    	fpop <<dt*record_time<< " " << pop_hist[dt] <<endl;
        fmark<<dt*record_time<< " "<<mark_hist[dt]<<endl;
    	//fhet <<t*record_time<< " " << het_hist[t] <<endl;

    }


    ///final out put
    
    time_t time_end;
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_), Fixed deme origin" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< run_time <<", "<<fixed_marker <<endl;

    cout << "Finished in " << t << " seconds \n";

    


    return 0;
}



