
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=450



g++ -o outf phage_inf_sim_het.cpp -lgsl
for b in 5 10 20 50 100 200
do
   for t in 50 100 200 500 1000 2000 5000
   do
      for (( i = 0; i < 10 ; i++ )) ### Outer for loop ###
      do
      	./outf -b $b -t $t -i $i &
      done 	
   done
done



