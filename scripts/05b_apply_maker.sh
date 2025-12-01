#!/bin/bash

#SBATCH --time=7-0
#SBATCH --mem=120G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=Maker
#SBATCH --output=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/Maker_gene_annotation_%j.out
#SBATCH --error=/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/log/Maker_gene_annotation_%j.err

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation"

REPEATMASKER_DIR="/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"
export PATH=$PATH:"/data/courses/assembly-annotation-course/CDS_annotation/softwares/RepeatMasker"

# ========================================
# Load modules
# ========================================
module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

# ========================================
# Run MAKER with MPI
# ========================================
echo "[INFO] Starting MAKER annotation..."
echo "[INFO] Start time: $(date)"

cd "$WORKDIR/results/02_gene_annotations"

mpiexec --oversubscribe -n 50 apptainer exec \
    --bind $SCRATCH:/TMP --bind $COURSEDIR --bind $REPEATMASKER_DIR \
    --bind /data \
    $COURSEDIR/containers/MAKER_3.01.03.sif \
    maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_evm.ctl maker_exe.ctl

echo "[INFO] End time: $(date)"