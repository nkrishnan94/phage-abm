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
const int N_demes = 150; // number of demes in comiving frame
//const int N_spec = 2; // number of 'species' including phage and bacteria
const int K_bac=100; // deme size for bacteria
const int K_vir = 100; // deme size for phage - >beta*K_bac*2
float tao = 200000; // lysis time in simulation steps
int beta = 30; //number of phage released with lysis
float M = 1; // Migration rate
int prof_hist = 0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 1*pow(10,6); // Run time in generations
unsigned int assign_gene_time = 5*pow(10,6); // Run time in generations
int samp_id=0;
float alpha = 1;
bool from_file_flag = false;



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


    string termination = "_demes" + strDeme.str() +"_" +strTime.str() + ".txt";
    string velName = "velocity_Nb" + strKb.str()  + "_migr" + strM.str() +"_B" +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str() +"_ID"+strI.str()+  termination;
    string hetName = "hetero_Nb" + strKb.str()  + "_migr" + strM.str() +"_B" +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str() +"_ID"+strI.str()+  termination;
    string profpName = "profile_phage_Nb" + strKb.str()  + "_migr" + strM.str()  +"_B" +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str()+"_ID"+strI.str()+ termination;
    string profbName = "profile_bac_Nb" + strKb.str()  + "_migr" + strM.str() +"_B" +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str() +"_ID"+strI.str()+ termination;
    string logName = "log_Nb" + strKb.str()  + "_migr" + strM.str() + "_B"  +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str()+"_ID"+strI.str()+  termination;
    //string markName = "markers_" + strKb.str()  + "_migr" + strM.str() + "_B" +strB.str()+"+_tau"+strT.str()+"_alpha"+strA.str()+ "_ID"+strI.str()+termination;
    string folder = "anc_data/";
    flog.open(folder+logName);
    //fhet.open(hetName);
    fpop.open(folder+velName);
    fprofp.open(folder+profpName);
    fprofb.open(folder+profbName);
    //fmark.open(folder+markName);


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
    //long V_tot[N_demes] = {0};
    //long V_prof[N_demes] = {0};
    double shiftDemes = 0; // Number of demes shifted
    int record_time=100000;
    vector <double> pop_hist;
    //vector <double> mark_hist;
    vector <double> het_hist;
    int total_phage = int(N_demes/2)*100;
    int total_bac_no_inf=0;
    int total_bac_inf= K_bac*int(N_demes/2);
    //int tao = tao_*K_bac*beta;
    long tao_count = pow(10,int(log10(tao) + 2));
    int t=0;
    int demes_alive=0;
    int fixed_marker;

    int shiftpop=0;
    int N_mark = 2;

    if(from_file_flag==true){
        string line;
        ifstream myfile ("fisher.txt");
        int j = 0;
        if (myfile.is_open())
        {


            while ( getline (myfile,line) )
            {
              string::iterator it;
              int index = 0;
              //for ( it = line.begin() ; it < line.end(); it++ ,index++)
              //{
               // cout << *it;
                //cout << line << '\n';
              //}
              V_deme[j][0] = stof(line);
              //cout << stoi(line);
              V_deme[j][1] = 0;
              j+=1;
              //cout << line[1,2,3];

           
            }
            myfile.close();


        }


    } else{

        for(int m= 0; m<int(N_demes/2);m++){

            V_deme[m][0]=100;


        }


    }

    ///---setup iinitial population
    for(int m= 0; m<int(N_demes/2);m++){

		V_deme[m][0]=100;

		//V_deme[m][1]=10;

    }


    ///-------Helper Funcs-----------


    //main loop
    //while((t<1.1*assign_gene_time) ){
    while((demes_alive!=1) ||(t<assign_gene_time) ){
    //for (unsigned int t = 0; t < N_gen; t++){
        if (total_phage>1){
            

            //migrationt
            uniform_int_distribution<int> distribution_ind(0, total_phage-1);
            int phage_ind = distribution_ind(e)+1;
            //uniform_int_distribution<int> distribution_Im(0, K_vir);
            //int r_ind = distribution_Im(e);
            int phage_cnt=0;
            int phage_found=0;
            int r_deme;
            int r_phage;
            for(int m =0;m<N_demes;m++){
                for(int n =0;n<N_mark;n++){

                    if (phage_found==0){

                        phage_cnt+=V_deme[m][n];

                        if (phage_cnt+1>phage_ind){
                            r_deme=m;
                            if(t>assign_gene_time){
                                r_phage = n;


                            } else{

                                r_phage = 0;

                            }
            
                            phage_found=1;

                        }
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
            
            uniform_real_distribution<double> distribution_d(0.0, 1.0);
            double p_mig = distribution_d(e);
            //cout<<r_deme<<endl;


            if(r_deme>0){
                if(p_mig<M){
                    double r_dir = distribution_d(e);
                    if (r_dir<.5){
                        mig_deme = r_deme -1;
                    } else{
                        mig_deme = r_deme +1;

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
            uniform_int_distribution<int> distribution_indi(0, total_phage-1);
            int phage_indi = distribution_indi(e) +1;
            uniform_real_distribution<double> distribution_di(0.0, 1.0);
            double p_inf = distribution_di(e);
            int phage_cnti=0;
            int phage_foundi=0;
            int r_demei;
            int r_phagei;
            if (p_inf<alpha){
                for(int m =0;m<N_demes;m++){
                    for(int n =0;n<N_mark;n++){

                        if (phage_foundi==0){

                            phage_cnti+=V_deme[m][n];
                            if (phage_cnti+1>phage_indi){
                                r_demei=m;

                                if(t>assign_gene_time){
                                    r_phagei = n;


                                } else{

                                    r_phagei = 0;

                                }
 
                                phage_foundi=1;

                            }

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
        //cout<<tot_pop<<endl;
        //cout <<tot_pop<<endl;

        ///segmentation fault here!!!
        int shift = 0;
        //if (last_deme >= N_demes*(3/4)){
        shift = last_deme-5-int(N_demes/2);
            
        
        //cout<<last_deme<<endl;
        

        if (shift>0){

            //cout<<shift<<endl;
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
            for(int m=0;m<N_demes;m++){
                //int deme_tot =0;
                int org_tot=0;

                for(int n=0;n<N_mark;n++){


                    org_tot+=V_deme[m][n];
                    //deme_tot+=V_deme[m][n];

                }
                cout<< org_tot <<" ";
                //cout<< deme_tot <<" ";
                        //for(int n=0;n<N_demes;n++){
        //    V_tot[n]=0;
        //}



            }

            cout<<endl;

            for(int m=0;m<N_demes;m++){
                int deme_tot =0;
                //int org_tot=0;

                for(int n=0;n<N_mark;n++){


                    //org_tot+=V_deme[n][m];
                    deme_tot+=V_deme[n][m];

                }
                //cout<< org_tot <<" ";
                cout<< deme_tot <<" ";
                        //for(int n=0;n<N_demes;n++){
        //    V_tot[n]=0;
        //}



            }
            cout<<endl;




            //het = calcHet(V_deme);
            //cout<<"timstep: "<< t<<" Het: "<< het<<endl;

            cout<<"timstep: "<< t <<" # Markers alive: "<< demes_alive<<endl;
            cout<<"Demes until shift: "<< shift<<endl;
            cout<<"Last demes: "<< last_deme<<endl;
            cout<<"shift pop: "<< shiftpop <<" in bounds pop "<<total_phage<<endl;
            pop_hist.push_back(shiftpop+total_phage);
            //mark_hist.push_back(demes_alive);
            

            //het_hist.push_back(het);



            
            /*if (t>assign_gene_time){

                //cout<<"Saving profiles\n";

                //ostringstream strDT;

                //strDT << t;

                //for(int m=0;m<N_demes;m++){
                 //   cout<<V_tot[m]<<" ";



                //}
                //cout <<endl;

                string proftpName = "profp_T_"+strDT.str()+" " + strTime.str() + ".txt";
                string proftbName = "profb_T_"+strDT.str()+" " + strTime.str() + ".txt";
                ofstream fproftp,fproftb;
                fproftp.open(folder+proftpName);
                fproftb.open(folder+proftbName);
                for (int m = 0; m < N_demes; m++){
                    fproftp << m <<" ";

                    for (int n=0;n<N_mark;n++){

                        fproftp <<" "<<V_deme[m][n];


                    }
                    fproftp <<endl;

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
                    fproftb <<m<<" "<< B_health <<" " <<B_inf<<" " <<B_lys <<endl;


                }



            }*/






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
    
    //std::cout <<"hi"<<std::endl;
    for(int dt=0; dt <int(t/record_time);dt++){

    	fpop <<dt*record_time<< " " << pop_hist[dt] <<endl;
        //fmark<<dt*record_time<< " "<<mark_hist[dt]<<endl;

    	//fhet <<t*record_time<< " " << het_hist[t] <<endl;

    }


    ///final out put
    
    time_t time_end;
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_), Fixed deme origin" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< run_time <<", "<<fixed_marker <<endl;

    cout << "Finished in " << run_time << " seconds \n";

    


    return 0;
}



