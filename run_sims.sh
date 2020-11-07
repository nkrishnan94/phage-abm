
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=450



g++ -o outf phage_inf_sim_het.cpp -lgsl
for b in 10 20 50 100 200
do
   for t in 5000 10000 20000 50000 100000
   do
      for (( i = 0; i < 1 ; i++ )) ### Outer for loop ###
      do
      	./outf -b $b -t $t &
      done 	
   done
done



