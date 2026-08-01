#!/bin/bash
#SBATCH --job-name=dvi-paper-sim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --requeue

# Phần nặng nhất của profile `paper`: mô phỏng chính (Kịch bản 1-5, 3 nuisance
# learner, 100 lần lặp) và benchmark baseline. Hai stage này phải chạy cùng nhau
# vì bảng so sánh cần kết quả psi_0 của mô phỏng chính.
# Nộp lại nhiều lần cho tới khi log báo đủ số tổ hợp; mỗi lần chạy tiếp phần còn
# thiếu nhờ checkpoint trong results/R/ và checkpoints/R/.

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs
exec > "logs/paper-sim-${SLURM_JOB_ID}.out" 2>&1
module purge 2>/dev/null || true
module load R 2>/dev/null || true
source ~/.bashrc

echo "Job $SLURM_JOB_ID on $(hostname); started $(date)"

Rscript main.R \
    --profile paper \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --project-root "$SLURM_SUBMIT_DIR" \
    --stages simulations,baseline

echo "Finished $(date)"
