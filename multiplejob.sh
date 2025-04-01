#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --time=3-00
#SBATCH --mem=244G
#SBATCH --array=0-1
#SBATCH --qos=bbdefault
#SBATCH --output=./slurm_logs/slurm-%j.out
#SBATCH --mail-type=ALL


FILENAME_LIST=($(<i.txt))  # Creates an indexed array from the contents of input_list.txt
INPUT_FILENAME=${FILENAME_LIST[${SLURM_ARRAY_TASK_ID}]}  # Look-up using array index

echo "I am array index ${SLURM_ARRAY_TASK_ID} and am processing file: ${INPUT_FILENAME}"

nextflow run main_hazex_v.2.nf --paired_reads="${INPUT_FILENAME}/*{1,2}.fq.gz" \
--reference_genome="/rds/projects/h/hejq-msc-2024-wgbs/Test_Sbatch/data/references/Hazelnut_CavTom2PMs-1.0/GCF_932294415.1" \
--reference_name="GCF_932294415.1_dhQueRobu3.1_genomic.fa" --index_requirement=0 \
--pipeline_loc="/rds/projects/h/hejq-msc-2024-wgbs/Test_Sbatch"