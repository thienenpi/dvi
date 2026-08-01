#!/bin/bash
#SBATCH --job-name=dvi-smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00

# Chạy nhanh để kiểm tra môi trường, đường dẫn dữ liệu và toàn bộ đường đi của
# pipeline (n = 500, 2 lần lặp, 2 kịch bản). Kết quả không dùng cho báo cáo.

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs
exec > "logs/smoke-${SLURM_JOB_ID}.out" 2>&1
module purge 2>/dev/null || true
source ~/.bashrc
conda activate dvi

echo "Job $SLURM_JOB_ID on $(hostname); started $(date)"
Rscript -e 'cat(R.version.string, "\n")'

Rscript main.R \
    --profile smoke \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --project-root "$SLURM_SUBMIT_DIR"

echo "Finished $(date)"
