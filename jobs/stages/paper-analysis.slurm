#!/bin/bash
#SBATCH --job-name=dvi-paper-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --requeue

# Phần còn lại của profile `paper`: giá trị mục tiêu giải tích, ablation, hai bộ
# dữ liệu thực, metadata môi trường và bảng đối chiếu yêu cầu đồ án.
# Chạy sau khi jobs/stages/paper-simulations.slurm đã hoàn tất.
#
# Lưu ý: stage `check` đọc lại kết quả từ `state` trong cùng tiến trình, nên khi
# tách giai đoạn như thế này, các mục phụ thuộc mô phỏng chính trong bảng đối
# chiếu sẽ ở trạng thái "Chưa đạt". Bảng đầy đủ lấy từ lần chạy đủ stage
# (jobs/reproduce/paper.slurm), còn ở đây dùng để kiểm tra phần dữ liệu thực.

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs
exec > "logs/paper-analysis-${SLURM_JOB_ID}.out" 2>&1
module purge 2>/dev/null || true
source ~/.bashrc
conda activate dvi
Rscript scripts/install_packages.R

echo "Job $SLURM_JOB_ID on $(hostname); started $(date)"

Rscript main.R \
    --profile paper \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --project-root "$SLURM_SUBMIT_DIR" \
    --stages targets,ablation,real,report

echo "Finished $(date)"
