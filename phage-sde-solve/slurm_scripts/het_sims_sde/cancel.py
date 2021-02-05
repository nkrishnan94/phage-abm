import os
import glob
fs =glob.glob('*.log')[0]
f = open(fs,'r')
print('file: '+fs)
lines = f.readlines()

for line in lines:
	try:
	    num = line.split('job ')[1].split(' ')[0]
	    os.system('scancel '+str(num))
	    os.system('qstat '+str(num))
	except:
		None
