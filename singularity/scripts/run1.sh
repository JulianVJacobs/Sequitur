#!/bin/bash

# List of job names, node lists, and script inputs
declare -A jobs
jobs=(
    ["sequitur.50.0.2585"]="mscluster59 sequitur.py 50.0.2585"
    ["sequitur.80.0.1616"]="mscluster59 sequitur.py 80.0.1616"
    ["sequitur.100.0.1293"]="mscluster59 sequitur.py 100.0.1293"
    ["sequitur.300.0.431"]="mscluster59 sequitur.py 300.0.431"
    ["sequitur.50.1.2585"]="mscluster52 sequitur.py 50.1.2585"
    ["sequitur.80.1.1616"]="mscluster52 sequitur.py 80.1.1616"
    ["sequitur.100.1.1293"]="mscluster52 sequitur.py 100.1.1293"
    ["sequitur.300.1.431"]="mscluster52 sequitur.py 300.1.431"
    ["dbg.50.0.2585"]="mscluster53 dbg.py 50.0.2585"
    ["dbg.80.0.1616"]="mscluster53 dbg.py 80.0.1616"
    ["dbg.100.0.1293"]="mscluster53 dbg.py 100.0.1293"
    ["dbg.300.0.431"]="mscluster53 dbg.py 300.0.431"
    ["dbg.50.1.2585"]="mscluster54 dbg.py 50.1.2585"
    ["dbg.80.1.1616"]="mscluster54 dbg.py 80.1.1616"
    ["dbg.100.1.1293"]="mscluster54 dbg.py 100.1.1293"
    ["dbg.300.1.431"]="mscluster54 dbg.py 300.1.431"
)

# Path to the SLURM script template
SLURM_SCRIPT_TEMPLATE="job.slurm"

# Iterate through the jobs and submit each one
for job_name in "${!jobs[@]}"; do
    # Extract node list and script input
    IFS=' ' read -r node_list script_name inputs <<< "${jobs[$job_name]}"

    # Path to the modified SLURM script
    MODIFIED_SLURM_SCRIPT="modified_${job_name}.slurm"

    # Copy the template to a new file
    cp $SLURM_SCRIPT_TEMPLATE $MODIFIED_SLURM_SCRIPT

    # Use sed to replace the job name, node list, and script input in the SLURM script
    sed -i "s/^#SBATCH -J .*/#SBATCH -J $job_name/" $MODIFIED_SLURM_SCRIPT
    sed -i "s/^#SBATCH --nodelist=.*/#SBATCH --nodelist=$node_list/" $MODIFIED_SLURM_SCRIPT

    # Submit the modified SLURM script
    sbatch $MODIFIED_SLURM_SCRIPT $script_name $inputs
done