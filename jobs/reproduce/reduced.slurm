#!/bin/bash
#SBATCH --job-name=dvi-reduced
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --requeue

# Profile giảm tải: 10 lần lặp mô phỏng chính, 10 lần lặp baseline, 30 lần lặp
# ablation. Dùng để kiểm tra toàn bộ pipeline trước khi chạy profile `paper`.
# Kết quả được ghi tăng dần vào results/R/simulation_reduced.csv nên nộp lại
# đúng script này sẽ chạy tiếp từ tổ hợp còn thiếu thay vì làm lại từ đầu.

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs
exec > "logs/reduced-${SLURM_JOB_ID}.out" 2>&1
module purge 2>/dev/null || true
source ~/.bashrc
conda activate dvi

echo "Job $SLURM_JOB_ID on $(hostname); started $(date)"
Rscript -e 'cat(R.version.string, "\n")'

Rscript main.R \
    --profile reduced \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --project-root "$SLURM_SUBMIT_DIR"

echo "Finished $(date)"
