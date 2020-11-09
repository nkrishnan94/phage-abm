
#!/bin/bash
# A sample Bash script, by Ryan
cnt=0
tot=450



g++ -o outf50 phage_inf_sim_het.cpp -std=c++11 -lgsl

for t in .1 .2 .5 5 10
do
  for (( i = 0; i < 1 ; i++ )) ### Outer for loop ###
  do
    ./outf50 -t $t -i $i
  done
done




