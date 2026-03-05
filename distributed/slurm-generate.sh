#!/bin/bash
# This is a bash script which sets up a slurm job submission, submits it and also sets necessary enviroments for running exciting with multiple cores

chunk=$1

workdir="$HOME/opt/Phunky.jl"
cd "$workdir" || exit

#Create slurm log directory
mkdir -p "$workdir/distributed/job-submission"

#==================================================================================================================================
#==================================================CONTENTS OF THE run.slurm FILE==================================================
#==================================================================================================================================
cat > run.slurm << EOF
#!/bin/bash

#SBATCH -vvv
#SBATCH --job-name=calc_$chunk
#SBATCH --partition=short

# ==================CONTROL ALLOCATION==================
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
# ==================CONTROL ALLOCATION==================

#SBATCH --time=48:00:00             # Time limit hrs:min:sec
#SBATCH --output=${workdir}/distributed/job-submission/slurmjob_calc_$chunk.log    

echo "[SLURM-SETUP] Beginning to load necessary modules..."
julia --project=. --version

julia --project=. --threads=20 $workdir/distributed/calc_$chunk.jl

echo " "
echo "[SLURM-SETUP] End of slurm-script instructions!"
date

EOF
#==================================================================================================================================
#==================================================CONTENTS OF THE run.slurm FILE==================================================
#==================================================================================================================================

chmod 750 run.slurm
sbatch run.slurm

rm run.slurm

cd ${workdir}/distributed
