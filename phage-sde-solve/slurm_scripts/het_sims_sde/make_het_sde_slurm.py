import numpy as np
import glob
import os 
from datetime import datetime

alphas = np.array([.005,.01,.03])[::-1] #alpha values 
Bs = np.array([25,50,75,100])#B0 values (need corresponding script for each value in dir. )
tausl =4
taus=np.zeros((len(alphas),tausl))
aBt_0 = 50*.03
aBt_f = 200*.03
for a,alpha in enumerate(alphas):
    taus[a] = np.linspace((aBt_0)/alpha,(aBt_f)/alpha,tausl)
taus=taus.astype(int)
samps = 500
files= []
time_str=datetime.now()
flogstr ="submit_log"+".log"
flog = open(flogstr,"a")
flog.write('date/time:\n')
flog.write(str(time_str)+'\n')
flog.write('taus:\n')
flog.write(str(taus)+'\n')
flog.write('alphas:\n')
flog.write(str('alphas')+'\n')
flog.write('Bacs.:\n')
flog.write(str('Bs')+'\n')


for a,alpha in enumerate(alphas):
	for tau in taus[a]:
		for B in Bs:
			param_str="tau"+str(tau)+"_a"+str(alpha)+"_Nb"+str(B)
			f = open("run_het_sims_sde_"+param_str+".slurm","a")
			files.append("run_het_sims_sde_"+param_str+".slurm")
			f.write('#!/bin/bash\n')
			f.write('#SBATCH -J run_sims_Nb100'+param_str+'\n')
			f.write('#SBATCH -A FUSCO-SL3-CPU\n')
			f.write('#SBATCH --nodes=1\n')
			f.write('#SBATCH --ntasks=1\n')
			f.write('#SBATCH --cpus-per-task=1\n')

			if tau <301:
				f.write('#SBATCH --time=00:02:00\n')
			elif tau <601:
					f.write('#SBATCH --time=00:06:00\n')

			else:
				f.write('#SBATCH --time=00:12:00\n')


			f.write('#SBATCH --mem=5980mb\n')
			f.write('#SBATCH --array=1-'+str(samps)+'\n')
			f.write('#SBATCH -p skylake\n')
			f.write('./etc/profile.d/modules.sh\n')
			f.write('module purge \n')
			f.write('module load rhel7/default-peta4\n')
			f.write('echo "This is job" $SLURM_ARRAY_TASK_ID\n')
			f.write('g++ -o outanc_'+param_str+'_$SLURM_ARRAY_TASK_ID phage_inf_coarse_het_Nb'+str(B)+'.cpp -std=c++11\n')
			f.write('./outanc_'+param_str+'_$SLURM_ARRAY_TASK_ID -t '+str(tau)+' -a '+str(alpha)+' -i $SLURM_ARRAY_TASK_ID\n')
			f.close()

for f in files:
	tau = f.split('tau')[1].split('_')[0]
	alpha = f.split('_a')[1].split('_')[0]
	B = f.split('_Nb')[1].split('.')[0]
	#print(alpha)
	flog.write("Tau: "+str(tau) +" Alpha: " +str(alpha)+ " Bac. in deme: " +str(B)+'\n')
	os.system('sbatch '+str(f)+ ' >> '+flogstr)
os.system('mkdir slurm_scripts\n')
os.system('mv *.slurm slurm_scripts/\n')

		


