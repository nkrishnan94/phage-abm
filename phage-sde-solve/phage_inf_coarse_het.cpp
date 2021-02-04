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
#include <algorithm>


//define key parameters
const int N_demes = 200; // number of demes in comiving frame
//const int N_spec = 2; // number of 'species' including phage and bacteria
const int K_bac=100; // deme size for bacteria
const int K_vir = 100; // deme size for phage - >beta*K_bac*2
float tao = 800; // lysis time in simulation steps
int beta = 50; //number of phage released with lysis
float M = .25; // Migration rate
int prof_hist =  0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 1*pow(10,4); // Run time in generations
int samp_id=0;
float alpha = 0.03;
unsigned int assign_gene_time = 5*pow(10,4); // Run time in generations
//int mig_det_flag = 0;  // =1 binomial sampling done to ddetermine proportion from each allele migaratoin; =0, done determinisitally for total and each allele. 
//int round_flag = 1; // =1, migration proportions are simply rounded down, otherwise binomieal sampling done to determine whether rounding occurs




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
        else if (c=='a')
            alpha = atof(optarg);
        else if (c == 't')
            tao = atof(optarg); // migration probability
        else if (c == 'i')
            samp_id = atoi(optarg); //keep track of profile through time
        else if (c == 'H')
            prof_hist = atoi(optarg); //keep track of profile through time 
    }
    
    //----------initialize time-------------
    time_t time_start;
    clock_t c_start = clock(); // Stop clock
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
    string folder = "het_data_sde/";
    //string folder = "";
    flog.open(folder+logName);
    fhet.open(folder+hetName);
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
    mt19937 e(samp_id);

    ///------additional parameters
    long B_deme[N_demes][K_bac] = {{0}};// Keep track of Bacteria -> Healthy, infected, lysed
    long V_deme[N_demes][2] = {{0}};// //Keep track of N_spec species of phage
    long V_deme_aux[N_demes][2] = {{0}};
    double shiftDemes = 0; // Number of demes shifted
    int shiftpop=0;
    int record_time_het=1;
    int record_time_pop=5000;
    vector <double> pop_hist;
    vector <double> het_hist;
    int total_phage = int(N_demes/2)*100;
    int total_bac_no_inf=0;
    int total_bac_inf= K_bac*int(N_demes/2);
    int t=0;
    float avgH;
    //int tao = tao_*K_bac*beta;
    long tao_count = pow(10,int(log10(tao) + 2));
    uniform_real_distribution<double> distribution_d(0.0, 1.0);



    ///---setup iinitial population
    for(int m= 0; m<int(N_demes/2);m++){

        V_deme[m][0]=100;



        V_deme[m][1]=0;

    }


    ///-------Helper Funcs-----------


    //main loop
    int aux_p0;
    int aux_p1;
    int hold0;
    int hold1;

    //for (unsigned int t = 0; t < N_gen; t++){
    while((avgH>.0001)||(t<1.1*assign_gene_time) ){

        //cout<<"hi"<<endl;
        //cout<< 199<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;
        int m=0;



        //first demem
        aux_p0 =V_deme[m][0]; 
        aux_p1 =V_deme[m][1]; 
        //cout<<V_deme[0][1]<<endl;

        //float M0 =((1-M/2)* V_deme[m][0] + (M/2)*V_deme[m+1][0]);
        //float M1 = ((1-M/2)* V_deme[m][1] + (M/2)*V_deme[m+1][1]);
        float M0 =((1-M/2)* V_deme[m][0] + (M/2)*V_deme[m+1][0]);
        float M1 = ((1-M/2)* V_deme[m][1] + (M/2)*V_deme[m+1][1]);
        int Mtot = int(M0+M1);
        //cout<<M0/(M0+M1) <<endl;
        binomial_distribution<int> distribution_M0( Mtot,  M0/(M0+M1) );
        //binomial_distribution<int> distribution_M1( 1,float(M1 - int(M1) ));
        V_deme[m][0] = distribution_M0(e);
        V_deme[m][1] = Mtot - V_deme[m][0];

        //V_deme[m][1] = int(M1) + distribution_M1(e);

    


        //absorption

        //cout<< 199<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;
        int Bempty =0;
        for (int n=0;n<K_bac;n++){
            if (B_deme[m][n]==0){
                //if(B_deme[m][n]==0){
                Bempty+=1;
                //}
 
            }

        }
        //cout<< Bempty<<endl;

        if (((V_deme[m][0]+V_deme[m][1] )>0)&& (Bempty>0)){
            //cout<< m<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;
            
            

            binomial_distribution<int> distribution_0( (V_deme[m][0] +V_deme[m][1])*Bempty ,alpha/ K_bac);
            

            int atot = distribution_0(e);

            atot=min(atot,Bempty);
            //cout<< Bempty<<endl;




            binomial_distribution<int> distribution_1( atot,float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) );
            int a0 = distribution_1(e);
            int a1 = atot-a0;
            //a1  >0 when allele1 is empty
            //cout<< V_deme[m][0]<<" "<<V_deme[m][1]<<" "<<V_deme[m][0] / (float( V_deme[m][0])+float(V_deme[m][1])) <<endl;

            //cout<<M0/(M0+M1)<<" "<<a0<<" "<<a1<<" "<<atot<<endl;

           //cout<<a0<<endl;
            //cout<< 0<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;
            int Babs = 0;
            int n=0;
            while((Babs< a0 )&&(n<K_bac)){
                if (B_deme[m][n] ==0){
                    B_deme[m][n] = (1) * tao_count;
                    Babs+=1;
                    //cout<< n<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;

                }
                n+=1;


            }
            while((Babs< (a0+a1) )&&(n<K_bac)){
                if (B_deme[m][n] ==0){
                    B_deme[m][n] = (2) * tao_count;
                    Babs+=1;
                    //cout<< n<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;

                }
                n+=1;


            }



            if (Babs<(a0+a1)){


                cout<<"Timestep: "<< t << " Deme: "<< m<<"total bac. "<< K_bac<<" "<<"Infections done"<< Babs<<" "<<"Bempty: "<< Bempty<< " Bac. Index: " <<a0<<" +  "<<a1<<endl;
                cout<<
                cout <<"NO UNINFECTED BACTERIA AVAILABLE FOR ABSORPTION" <<t <<endl;
                for (int n=0;n<K_bac;n++){
                    //if (B_deme[m][n] ==0){
                    cout<<n<<" "<<B_deme[m][n] <<endl;

                }
                exit(EXIT_FAILURE);
                
            }
            //cout<< 199<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;
        }
        //cout<< 199<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;



        /*if (V_deme[m][0]+V_deme[m][1]>0){
            binomial_distribution<int> distribution_0( V_deme[m][0],alpha*float(V_deme[m][0]/ float(V_deme[m][0]+V_deme[m][1]) ) );
            binomial_distribution<int> distribution_1( V_deme[m][1],alpha*float(V_deme[m][1]/ float(V_deme[m][0]+V_deme[m][1]) ) );

            int a0 = distribution_0(e);
            int a1 = distribution_1(e);
            //cout<<a0<<endl;


            
            int bac_ind;

            vector <int> phage_abs0((a0+a1),0);
            if (phage_abs0.size()>0){

                for(int v0=0; v0<a0;v0++){
                    phage_abs0[v0]=0;

                }
                for(int v1=a0; v1<a0+a1;v1++){
                    phage_abs0[v1]=1;

                }
                shuffle(begin(phage_abs0), end(phage_abs0),default_random_engine(samp_id));
                int r_phagei;

                



                for(int v=0; v<(a0+a1);v++){
                    r_phagei = phage_abs0[v];

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
        }*/
        


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
        //cout<< 199<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;
    
        for(int m=1;m<N_demes-2;m++){
            
            hold0 = aux_p0;
            hold1 = aux_p1;
            aux_p0 =V_deme[m][0]; 
            aux_p1 =V_deme[m][1]; 
            float M0 = ((1-M)* V_deme[m][0] + (M/2)*V_deme[m+1][0] + (M/2)*hold0);
            float M1 = ((1-M)* V_deme[m][1] + (M/2)*V_deme[m+1][1] + (M/2)*hold1);
            int Mtot = int(M0+M1);
            binomial_distribution<int> distribution_M0( Mtot,  M0/(M0+M1) );
            //binomial_distribution<int> distribution_M1( 1,float(M1 - int(M1) ));

            V_deme[m][0] = distribution_M0(e);
            V_deme[m][1] = Mtot -  V_deme[m][0];

            //cout <<m<<" "<<V_deme[m][0]<<" "<<V_deme[m][1]<<endl;

            //float M0 = ((1-M)* V_deme[m][0] + (M/2)*V_deme[m+1][0] + (M/2)*hold0);
            //float M1 = ((1-M)* V_deme[m][1] + (M/2)*V_deme[m+1][1] + (M/2)*hold1);
            //binomial_distribution<int> distribution_M0( 1,float(M0 - int(M0) ));
            //binomial_distribution<int> distribution_M1( 1,float(M1 - int(M1) ));

            //V_deme[m][0] = int(M0) + distribution_M0(e);
            //V_deme[m][1] = int(M1) + distribution_M1(e);

            int Bempty =0;
            for (int n=0;n<K_bac;n++){
                if (B_deme[m][n]==0){
                    if(B_deme[m][n]==0){
                        Bempty+=1;
                    }

                }

            }
            //cout<<Bempty<<endl;


            if ((V_deme[m][0]+V_deme[m][1]>0)&& (Bempty>0)){

                binomial_distribution<int> distribution_0( (V_deme[m][0] +V_deme[m][1])*Bempty ,alpha/ K_bac);
                

                int atot = distribution_0(e);

                binomial_distribution<int> distribution_1( atot,float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) );
                atot=min(atot,Bempty);
                int a0 = distribution_1(e);
                int a1 = atot-a0;


                //cout<<a0<<endl;
                int Babs=0;
                int n=0;
                while((Babs< a0 )&&(n<K_bac)){
                    if (B_deme[m][n] ==0){
                        B_deme[m][n] = (1) * tao_count;
                        Babs+=1;
                        //cout<< n<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;

                    }
                    n+=1;


                }
                while((Babs< (a0+a1) )&&(n<K_bac)){
                    if (B_deme[m][n] ==0){
                        B_deme[m][n] = (2) * tao_count;
                        Babs+=1;
                        //cout<< n<<" "<< V_deme[199][0]<<" "<<V_deme[199][1]<<endl;

                    }
                    n+=1;


                }
                if (Babs<(atot)){


                    cout<< "Timestep: "<< t<<" Deme: "<< m<<"Bempty: "<< Bempty<< "a0 a1: " << a0 <<" "<<a1 <<endl;
                    cout << "V0 "<< V_deme[m][0]<<" Abs* V_0* B " << alpha*float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) * float(Bempty)/float(K_bac) <<endl;

                    cout << "V1 "<< V_deme[m][1]<<"Abs* V_1* B " << alpha*float(V_deme[m][1] / float(V_deme[m][0]+V_deme[m][1]) ) * float(Bempty)/float(K_bac) <<endl;
                    cout <<"NO UNINFECTED BACTERIA AVAILABLE FOR ABSORPTION" <<t <<endl;
                    exit(EXIT_FAILURE);
                    
                }
            }

            /*if (V_deme[m][0]+V_deme[m][1]>0){
                //cout<<V_deme[m][0]+V_deme[m][1]<<endl;
                //absorption
                binomial_distribution<int> distribution_0( V_deme[m][0],alpha*float(V_deme[m][0]/ float(V_deme[m][0]+V_deme[m][1]) ) );
                binomial_distribution<int> distribution_1( V_deme[m][1],alpha*float(V_deme[m][1]/ float(V_deme[m][0]+V_deme[m][1]) ) );
                int a0 = distribution_0(e);
                int a1 = distribution_1(e);
                //cout << alpha*(V_deme[m][0]/ (V_deme[m][0]+V_deme[m][1]) )<< endl;

                
                int bac_ind;

                
                vector <int> phage_abs((a0+a1),0);
                if (phage_abs.size()){
                    //cout <<"absorption"<<endl;
                    for(int v0=0; v0<a0;v0++){
                    phage_abs[v0]=0;

                    }
                    for(int v1=a0; v1<a0+a1;v1++){
                        phage_abs[v1]=1;

                    }
                    shuffle(begin(phage_abs), end(phage_abs),default_random_engine(samp_id));
                    int r_phagei;
                    if ((a0+a1)>0){
                        //cout<<(a0+a1)<<endl;

                    }

                    for(int v=0; v<(a0+a1);v++){
                        r_phagei = phage_abs[v];
                        //cout<< r_phagei<<endl;
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
            }*/
            


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
            //cout <<m<<" "<<V_deme[m][0]<<" "<<V_deme[m][1]<<endl;



        }

        
        

        //shift populatoion
        total_phage=0;
        int last_deme=0;
        for(int m=0;m<N_demes;m++){


            total_phage+=V_deme[m][0]+V_deme[m][1];
            //cout<< m<<" "<< V_deme[m][0]<<" "<<V_deme[m][1]<<endl;
            if ((V_deme[m][0]+V_deme[m][1])>0){
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
        //cout<<last_deme<<" "<<shift<<endl;
            
        
        //cout<<last_deme<<endl;


        if (shift>0){
            //cout<<last_deme<<" "<<shift<<endl;
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



        if(t%record_time_pop ==0){

            cout<<endl;


            for(int m=0;m<N_demes;m++){
                //int total_deme=0;


                cout<<V_deme[m][0] <<" ";

            
            }
            cout<<endl;

            for(int m=0;m<N_demes;m++){
                //int total_deme=0;


                cout<<V_deme[m][1] <<" ";

            
            }
            cout<<endl;

            //avgH =calcHet(V_deme);

            cout<<"timestep: "<< t<<" Het: "<<avgH<< "total pop "<< total_phage<<endl;
            cout<<"demes until shift"<< shift <<endl;
            pop_hist.push_back(shiftpop+total_phage);
            het_hist.push_back(calcHet(V_deme));
        }

        if((t%record_time_het ==0)&&(t>assign_gene_time)){

            avgH =calcHet(V_deme);
            cout << "Average Het. :" <<avgH<<endl;
            het_hist.push_back(avgH);



        }
        if(t==assign_gene_time){


            

            cout<<"Assigning labels\n";
            // assign ing labels to free phage

            for(int m=0;m<N_demes;m++){

                int deme_pop = V_deme[m][0]+V_deme[m][1];
                if (deme_pop>0){

                    binomial_distribution<int> distribution_bi(deme_pop,0.5);
                    int a0_pick = distribution_bi(e);
                    V_deme[m][0] = a0_pick;
                    V_deme[m][1] = deme_pop - a0_pick;

                }

            
            }
            cout<<endl;


            for(int m=0;m<N_demes;m++){
                //int total_deme=0;


                cout<<V_deme[m][0] <<" ";

            
            }
            cout<<endl;

            for(int m=0;m<N_demes;m++){
                //int total_deme=0;


                cout<<V_deme[m][1] <<" ";

            
            }
            cout<<endl;

            //reset bacteria
            for(int m=0; m< N_demes;m++){
                for(int nb=0; nb< K_bac;nb++){
                        //uniform_real_distribution<double> distribution_dr(0.0, 1.0);
                        float p_all = distribution_d(e);


                        B_deme[m][nb] = B_deme[m][nb]*int(round(p_all)+1);
                        //cout<< round(p_all)+1<<endl;
                    

                }

            }

            cout<<"Saving profiles\n";
            for (int m = 0; m < N_demes; m++){
                


                fprofp << m <<" " <<V_deme[m][0]<<" " <<V_deme[m][1] <<endl;

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

            // assign ing labels to absorbed phaage


        }
        t+=1;

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
    for(int dt=0; dt <int(t/record_time_pop);dt++){

        fpop <<dt*record_time_pop<< " " << pop_hist[dt] <<endl;
        

    }

    for(int dt=0; dt <het_hist.size();dt++){


        //fpop <<dt*record_time<< " " << pop_hist[dt] <<endl;
        fhet <<assign_gene_time+dt*record_time_het<< " " << het_hist[dt] <<endl;

    }
    cout <<het_hist.size();


    ///final out put
    
    time_t time_end;
    clock_t c_fin = clock(); // Stop clock
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< double(c_fin - c_start)/CLOCKS_PER_SEC <<endl;

    cout << "Finished in " << double(c_fin - c_start)/CLOCKS_PER_SEC << " seconds \n";

    


    return 0;
}



