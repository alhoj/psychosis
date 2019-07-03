#!/bin/sh

chmod 755 /scratch/nbe/psykoosi/scripts/jobs/*

# This is the part where we submit the jobs that we cooked

for j in $(ls -1 "/scratch/nbe/psykoosi/scripts/jobs/");do
sbatch "/scratch/nbe/psykoosi/scripts/jobs/"$j
sleep 0.01
done
echo "All jobs submitted!"

#rm slurm-*
