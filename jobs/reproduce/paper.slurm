#!/bin/bash
#SBATCH --job-name=dvi-paper
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --requeue

# Profile của bài báo: 100 lần lặp, n = 10.000, lưới tương quan 0 -> 0,95, rừng
# hồi quy dùng số cây mặc định của grf. Đây là cấu hình để đối chiếu Table 2.
# Chi phí lớn, chủ yếu do ước lượng psi_0 (KDE và tích phân Monte Carlo).
# Job resume được: nộp lại script này sau khi hết thời gian sẽ chạy tiếp.
# Nếu tường thời gian của cluster ngắn hơn, dùng jobs/stages/ để tách giai đoạn.

cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs
exec > "logs/paper-${SLURM_JOB_ID}.out" 2>&1
module load R 2>/dev/null || true
source ~/.bashrc
conda activate dvi

echo "Job $SLURM_JOB_ID on $(hostname); started $(date)"
Rscript -e 'cat(R.version.string, "\n")'

Rscript main.R \
    --profile paper \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    --project-root "$SLURM_SUBMIT_DIR"

echo "Finished $(date)"
