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
const int K_bac=75; // deme size for bacteria
const int K_vir = 100; // deme size for phage - >beta*K_bac*2
float tao = 200; // lysis time in simulation steps
int beta = 50; //number of phage released with lysis
float M = .25; // Migration rate
int prof_hist =  0; // flag to keep track of history profile history through time, off by default
unsigned int N_gen = 1*pow(10,4); // Run time in generations
int samp_id=0;
float alpha = 0.03;
unsigned int assign_gene_time = 5*pow(10,4); // Run time in generations


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
    while ((c = getopt (argc, argv, "b:t:a:i:H")) != -1)
    {

        if (c == 'b')
            beta = atof(optarg); // migration probability
        else if (c=='a') //absorption
            alpha = atof(optarg);
        else if (c == 't')
            tao = atof(optarg); // lysis time
        else if (c == 'i')
            samp_id = atoi(optarg); //sample_id, to pass to seed
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
    mt19937 e(samp_id);

    ///------additional parameters
    long B_deme[N_demes][K_bac] = {{0}};// Keep track of Bacteria -> Healthy, infected, lysed
    long V_deme[N_demes][2] = {{0}};// //Keep track of N_spec species of phage
    long V_deme_aux[N_demes][2] = {{0}};
    double shiftDemes = 0; // Number of demes shifted
    int shiftpop=0; //initialize population shift. variable
    int record_time_het=5; //recording interval for heteozygoisty after assigning
    int record_time_pop=5000;
    vector <double> pop_hist; //initalize vector for population
    vector <double> het_hist; //initialize vector for avg. heterozygoisty
    int total_phage = int(N_demes/2)*100; //initial value for total_phage
    int t=0; //initialized time variable
    float avgH;  //initialize variable for current
    long tao_count = pow(10,int(log10(tao) + 2)); //
    int aux_p0; /// save unmigrated allele 0 pop for next deme
    int aux_p1; //save unmigrated allele 1 pop for next deme
    int hold0; //holder for prev. deme unmigrated allele 0 pop
    int hold1; //holder for prev. deme unmigrated allele 0 pop
    float H_thresh = 0.0001;
    uniform_real_distribution<double> distribution_d(0.0, 1.0); // initialize uniform dist. for 0-1.



    ///---setup iinitial population-----------
    for(int m= 0; m<int(N_demes/2);m++){

        V_deme[m][0]=100;
        V_deme[m][1]=0;

    }



    //---------main loop-----------
    while((avgH>H_thresh)||(t<1.1*assign_gene_time) ){

        //first deme
        int m=0;

         //MIGRATION
        //save unmigrated deme populations for next deme
        aux_p0 =V_deme[m][0]; 
        aux_p1 =V_deme[m][1]; 

        //find expected value of each allele after migration
        float M0 = ((1-M/2)* V_deme[m][0] + (M/2)*V_deme[m+1][0]);
        float M1 = ((1-M/2)* V_deme[m][1] + (M/2)*V_deme[m+1][1]);
        int Mtot = int(M0+M1); //total number of individuals after migration

        binomial_distribution<int> distribution_M0( Mtot,  M0/(M0+M1) ); // binomial sampling of allele 0 from total pop. based on expected fraction
        
        V_deme[m][0] = distribution_M0(e); //assign allele 0 from sampel number
        V_deme[m][1] = Mtot - V_deme[m][0]; //assign rest to allele

        //ABSORPTION
        //check for empty bacteria
        int Bempty =0;
        for (int n=0;n<K_bac;n++){
            if (B_deme[m][n]==0){
                Bempty+=1; // check all bacteria in deme, and add 1 if it is empty (value set to 0)
            }

        }
        //absorption can only possibly happen if there are available bacteria and phage
        if (((V_deme[m][0]+V_deme[m][1] )>0)&& (Bempty>0)){

            binomial_distribution<int> distribution_0( (V_deme[m][0] +V_deme[m][1])*Bempty , alpha/K_bac); //sample <alpha Bi Bi> from binomial distrib.

            int atot = distribution_0(e);
            atot=min(atot,Bempty); //atot : total phage to be absorbed (must be less than or equal to  Bempty)
            //cout<< Bempty<<endl;

            binomial_distribution<int> distribution_1( atot,float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) ); //binomially sample each allele to be absorbed based on allele fraction
            int a0 = distribution_1(e); //number of allele 0  to be absorbed
            int a1 = atot-a0; //assign rest to allele 1

            int Babs = 0; //keep track of how many absorption events occur
            int n=0;  //keep track of how any bacteria have been iterated through 

            //while there are bacteria available and the there are allele 0 phage left to be absorbed 
            while((Babs< a0 )&&(n<K_bac)){
                if (B_deme[m][n] ==0){
                    B_deme[m][n] = (1) * tao_count;
                    Babs+=1;
                }
                n+=1;
            }

            //while there are bacteria available and the there are allele 0 phage left to be absorbed
            while((Babs< (a0+a1) )&&(n<K_bac)){
                if (B_deme[m][n] ==0){
                    B_deme[m][n] = (2) * tao_count;
                    Babs+=1;
                }
                n+=1;
            }


            //ERROR IF NOT ALL PHAGE HAVE BEEN ABSORBED THAT ARE SCHEDUELED TO 
            if (Babs<(a0+a1)){

                cout<<"Timestep: "<< t << " Deme: "<< m<<"total bac. "<< K_bac<<" "<<"Infections done"<< Babs<<" "<<"Bempty: "<< Bempty<< " Bac. Index: " <<a0<<" +  "<<a1<<endl;
                cout<<
                cout <<"NO UNINFECTED BACTERIA AVAILABLE FOR ABSORPTION" <<t <<endl;
                for (int n=0;n<K_bac;n++){
                    cout<<n<<" "<<B_deme[m][n] <<endl;
                }
                exit(EXIT_FAILURE);
                
            }
        }


        //LYSIS
        for(int nb=0; nb< K_bac;nb++){
            if (B_deme[m][nb]>0){
                B_deme[m][nb]+=1; //advance all infected phage 'timer' by 1
            }

            //check if timer = tau for any
            if ((B_deme[m][nb]% tao_count)==tao){
                int burst_phage=(B_deme[m][nb]-tao)/tao_count-1; //find allele of new phage
                B_deme[m][nb]=-1; //set bacteria to 'dead'
                V_deme[m][burst_phage]+=beta; //add beta phage with  approiate allele to deme
            }

            if ((B_deme[m][nb]% tao_count)>tao){
                cout <<"UNLYSED BACTERIA AFTER LYSIS TIME!"<<endl;
                exit(EXIT_FAILURE);
            }


        }
        
        //subsequent demes
        for(int m=1;m<N_demes-2;m++){
            

            hold0 = aux_p0; //holder for previous deme, allele 0 
            hold1 = aux_p1; //holder for previous deme, allele 1
            aux_p0 =V_deme[m][0]; //assign current deme premigrated pop. for next deme, allele 0 
            aux_p1 =V_deme[m][1]; //assign current deme premigrated pop. for next deme, allele 1

            //find expected value of each allele after migration
            float M0 = ((1-M)* V_deme[m][0] + (M/2)*V_deme[m+1][0] + (M/2)*hold0);
            float M1 = ((1-M)* V_deme[m][1] + (M/2)*V_deme[m+1][1] + (M/2)*hold1);
            int Mtot = int(M0+M1);//total number of individuals after migration

            binomial_distribution<int> distribution_M0( Mtot,  M0/(M0+M1) ); // binomial sampling of allele 0 from total pop. based on expected fraction
            
            V_deme[m][0] = distribution_M0(e); //assign allele 0 from sampel number
            V_deme[m][1] = Mtot - V_deme[m][0]; //assign rest to allele


            //ABSORPTION
            //check for empty bacteria
            int Bempty =0;
            for (int n=0;n<K_bac;n++){
                if (B_deme[m][n]==0){
                    Bempty+=1; // check all bacteria in deme, and add 1 if it is empty (value set to 0)
                }

            }
            //absorption can only possibly happen if there are available bacteria and phage
            if (((V_deme[m][0]+V_deme[m][1] )>0)&& (Bempty>0)){

                binomial_distribution<int> distribution_0( (V_deme[m][0] +V_deme[m][1])*Bempty , alpha/K_bac); //sample <alpha Bi Bi> from binomial distrib.

                int atot = distribution_0(e);
                atot=min(atot,Bempty); //atot : total phage to be absorbed (must be less than or equal to  Bempty)
                //cout<< Bempty<<endl;

                binomial_distribution<int> distribution_1( atot,float(V_deme[m][0] / float(V_deme[m][0]+V_deme[m][1]) ) ); //binomially sample each allele to be absorbed based on allele fraction
                int a0 = distribution_1(e); //number of allele 0  to be absorbed
                int a1 = atot-a0; //assign rest to allele 1

                int Babs = 0; //keep track of how many absorption events occur
                int n=0;  //keep track of how any bacteria have been iterated through 

                //while there are bacteria available and the there are allele 0 phage left to be absorbed 
                while((Babs< a0 )&&(n<K_bac)){
                    if (B_deme[m][n] ==0){
                        B_deme[m][n] = (1) * tao_count;
                        Babs+=1;
                    }
                    n+=1;
                }

                //while there are bacteria available and the there are allele 0 phage left to be absorbed
                while((Babs< (a0+a1) )&&(n<K_bac)){
                    if (B_deme[m][n] ==0){
                        B_deme[m][n] = (2) * tao_count;
                        Babs+=1;
                    }
                    n+=1;
                }


                //ERROR IF NOT ALL PHAGE HAVE BEEN ABSORBED THAT ARE SCHEDUELED TO 
                if (Babs<(a0+a1)){

                    cout<<"Timestep: "<< t << " Deme: "<< m<<"total bac. "<< K_bac<<" "<<"Infections done"<< Babs<<" "<<"Bempty: "<< Bempty<< " Bac. Index: " <<a0<<" +  "<<a1<<endl;
                    cout<<
                    cout <<"NO UNINFECTED BACTERIA AVAILABLE FOR ABSORPTION" <<t <<endl;
                    for (int n=0;n<K_bac;n++){
                        cout<<n<<" "<<B_deme[m][n] <<endl;
                    }
                    exit(EXIT_FAILURE);
                    
                }
            }


            //LYSIS
            for(int nb=0; nb< K_bac;nb++){
                if (B_deme[m][nb]>0){
                    B_deme[m][nb]+=1; //advance all infected phage 'timer' by 1
                }

                //check if timer = tau for any
                if ((B_deme[m][nb]% tao_count)==tao){
                    int burst_phage=(B_deme[m][nb]-tao)/tao_count-1; //find allele of new phage
                    B_deme[m][nb]=-1; //set bacteria to 'dead'
                    V_deme[m][burst_phage]+=beta; //add beta phage with  approiate allele to deme
                }
            }

        }/// Steps for all demes completed


        

        //Shift population
        total_phage=0; // set total number of phage
        int last_deme=0; //initialize to keep track of the furtherest deme with phage
        for(int m=0;m<N_demes;m++){
            total_phage+=V_deme[m][0]+V_deme[m][1];
            if ((V_deme[m][0]+V_deme[m][1])>0){
                last_deme=m; //deme with largest m where phage total pop >0 is the 'a'

            }

            //Error for negative populations
            if(( V_deme[m][0]<0 )|| (V_deme[m][1]<0)){

                cout <<"deme population negative, timestep " <<t <<endl;
                cout <<V_deme[m][0]<<" " << V_deme[m][1]<< " "<<m<<endl;
                exit(EXIT_FAILURE);
            }
        
        }


        int shift = 0;
        shift = last_deme-10-int(N_demes/2); //calculate how many demes to shift by, to keep the population 10 demes past 'halway' mark

        //if shift is needed:
        if (shift>0){
            for(int s =0;s<shift;s++){
                shiftpop+= V_deme[s][0]+V_deme[s][1]; //keep track of number phage that will be lost to shift
            }
            

            for(int s =shift; s<N_demes;s++){
                V_deme[s-shift][0]=V_deme[s][0]; //shift allele 0 
                V_deme[s-shift][1]=V_deme[s][1]; //shift allele 1
                for(int b=0;b<K_bac;b++){
                    B_deme[s-shift][b]=B_deme[s][b]; //shift bacteria
                }
            }


            for(int s = N_demes-shift; s<N_demes;s++){
                V_deme[s][0]=0; //set previously occupied phage number to 0 
                V_deme[s][1]=0;
                for(int b=0;b<K_bac;b++){
                    B_deme[s][b]=0; //reset bacteria to all healthy
                }
            }


            
        }



        if(t%record_time_pop ==0){
            //print out populations for allele 0 
            cout<<endl;
            for(int m=0;m<N_demes;m++){
                //int total_deme=0;


                cout<<V_deme[m][0] <<" ";

            
            }
            cout<<endl;
            //print out populations for allele 1
            for(int m=0;m<N_demes;m++){
                //int total_deme=0;


                cout<<V_deme[m][1] <<" ";

            
            }
            cout<<endl;



            cout<<"timestep: "<< t<<" Het: "<<avgH<< "total pop "<< total_phage<<endl;
            cout<<"demes until shift"<< shift <<endl;
            pop_hist.push_back(shiftpop+total_phage); //record pop in vector
        }


        if((t%record_time_het ==0)&&(t>assign_gene_time)){
            avgH =calcHet(V_deme);
            cout << "Average Het. :" <<avgH<<endl;
            het_hist.push_back(avgH);
        }

        ///at specified time, aassign labels
        if(t==assign_gene_time){

            cout<<"Assigning labels\n";
            // assign ing labels to free phage
            for(int m=0;m<N_demes;m++){

                int deme_pop = V_deme[m][0]+V_deme[m][1];
                if (deme_pop>0){
                    binomial_distribution<int> distribution_bi(deme_pop,0.5); //binomially sample deme population where <f> =f
                    int a0_pick = distribution_bi(e);
                    V_deme[m][0] = a0_pick; //assign to allele 0
                    V_deme[m][1] = deme_pop - a0_pick; //give rest to other allele
                }
            
            }


            //print out profiles after assigning  labels
            cout<<endl;
            for(int m=0;m<N_demes;m++){
                cout<<V_deme[m][0] <<" ";
            }
            cout<<endl;

            for(int m=0;m<N_demes;m++){
                cout<<V_deme[m][1] <<" ";
            }
            cout<<endl;

            //absorbed phage
            for(int m=0; m< N_demes;m++){
                for(int nb=0; nb< K_bac;nb++){
                    if (B_deme[m][nb]>0){ //find absorbed phage
                        float p_all = distribution_d(e); //draw number from 0 -1
                        //THIS NEXT LINE ASSUMES THAT POPULATIONS START MONOALLELIC!!!!! OTHERWISE NEED TO EXTRACT TAO AND 
                        int curr_tao = (B_deme[m][nb]% tao_count);
                        B_deme[m][nb] = tao_count*int(round(p_all)+1) + curr_tao; //assign either aalele 0 or allele 1 with 50% probability (rounding random number)
                    }
                }

            }
            //save profile to file
            cout<<"Saving profiles\n";
            for (int m = 0; m < N_demes; m++){

                fprofp << m <<" " <<V_deme[m][0]<<" " <<V_deme[m][1] <<endl; //saving phage profiles

                //added up bacteria profiles for each bacterial state
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
                fprofb <<m<<" "<< B_health <<" " <<B_inf<<" " <<B_lys <<endl; //save bacterial profiles
            }

        }
        t+=1; //advance time by 1

    }
    
    //////-----write data files-------------
    for(int dt=0; dt <pop_hist.size();dt++){

        fpop <<dt*record_time_pop<< " " << pop_hist[dt] <<endl; //save population ffrom vector file 

    }

    for(int dt=0; dt <het_hist.size();dt++){

        fhet <<assign_gene_time+dt*record_time_het<< " " << het_hist[dt] <<endl; //save population ffrom vector file 

    }



    ///final output
    time_t time_end;
    clock_t c_fin = clock(); // Stop clock
    double run_time = difftime(time(&time_start), time(&time_end));
    flog << "Number of generations, Migration rate, Number of demes, Start time, Elapsed run time (secs_" << endl;
    flog << N_gen << " "  << M << ", " << N_demes << time_start<< double(c_fin - c_start)/CLOCKS_PER_SEC <<endl;
    cout << "Finished in " << double(c_fin - c_start)/CLOCKS_PER_SEC << " seconds \n";

    return 0;
}



