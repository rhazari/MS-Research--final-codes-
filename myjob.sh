rm *.pbs
for i in $(seq 1 1 100)
do
file="job_$i.pbs"
output="Algo_100_Min_$i.txt"   
echo  "#!/bin/bash" >> $file
echo  "#PBS -N Algo_Min" >> $file
echo  "#PBS -l select=1:ncpus=1,walltime=10:00:00" >> $file
echo  "#PBS -j oe" >> $file
echo  "cd \$PBS_O_WORKDIR" >> $file
echo  "rm Algo_Min.*" >> $file
echo  "python Algo1_sim_v1.py $output" >> $file
qsub $file
rm $file
done
