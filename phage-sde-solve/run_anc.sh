tot = 300
cnt = 0
for t in 50 100 200 500
do
	for (( i = 1; i <501; i++ )) ### Outer for loop ###
	do
		g++ -o outf_$t_$i phage_inf_coarse_anc.cpp
		./outf_$t_$i -t $t -i $i 
		cnt=$[cnt+1]
		echo $cnt "out of" $tot "\n"
	done 
done