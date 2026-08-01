# Tái thực nghiệm Decorrelated Variable Importance

Tái thực nghiệm phương pháp trong Verdinelli và Wasserman, *Decorrelated Variable
Importance* (JMLR 25, 2024 — bản PDF trong `docs/`): năm kịch bản mô phỏng của
bài báo, ba nuisance learner, thủ tục suy luận t-Cross, hai baseline đối chứng
(residual CPI, marginal PFI) và phân tích hai bộ dữ liệu thực Energy Efficiency
và Concrete Compressive Strength.

Hai bộ dữ liệu đã nằm sẵn trong `data/`, không cần tải thêm.

---

# Chạy trên máy cá nhân

Chạy lần lượt các lệnh dưới đây từ trên xuống.

## Bước 1 — Lấy mã nguồn

```bash
git clone https://github.com/thienenpi/dvi.git
cd dvi
```

## Bước 2 — Cài R

Bỏ qua bước này nếu `Rscript --version` đã chạy được. Chọn **một** cách:

```bash
# Ubuntu / Debian
sudo apt-get update && sudo apt-get install -y r-base build-essential

# hoặc conda (không cần quyền sudo)
conda create -y -n dvi -c conda-forge r-base && conda activate dvi

# hoặc macOS
brew install r
```

## Bước 3 — Cài package R

```bash
Rscript scripts/install_packages.R
```

Script cài `jsonlite`, `readxl`, `mgcv`, `grf` và in phiên bản từng package. Thêm
`--notebook` nếu muốn chạy notebook (`IRdisplay`, `IRkernel`).

## Bước 4 — Chạy thử toàn bộ pipeline

```bash
Rscript main.R --profile smoke --threads 4
```

Vài phút. Đây là bài kiểm tra môi trường, đường dẫn dữ liệu và toàn bộ đường đi
của pipeline; **kết quả không dùng cho báo cáo**.

## Bước 5 — Chạy cấu hình giảm tải

```bash
Rscript main.R --profile reduced --threads 8
```

10 lần lặp mô phỏng chính, 10 lần lặp baseline, 30 lần lặp ablation. Đây là cấu
hình đủ dùng cho phần lớn nội dung báo cáo.

## Bước 6 — Chạy cấu hình của bài báo (tùy chọn, rất tốn thời gian)

```bash
Rscript main.R --profile paper --threads 16
```

100 lần lặp, n = 10.000, lưới tương quan 0 → 0,95. Chỉ cấu hình này mới so sánh
trực tiếp được với Table 2 của bài báo. Nếu bị ngắt giữa chừng, chạy lại **đúng
lệnh trên** để tiếp tục phần còn thiếu.

## Bước 7 — Xem kết quả

```bash
ls results/R/    # bảng CSV
ls figures/R/    # hình PNG 300 DPI và PDF vector
cat README_EXPERIMENTS.md
```

---

# Chạy trên cluster Slurm

## Bước 1 — Đưa mã nguồn lên server

```bash
ssh <user>@<host> 'git clone https://github.com/thienenpi/dvi.git ~/dvi'
```

Hoặc đồng bộ từ máy cá nhân:

```bash
rsync -az --partial --info=progress2 \
    --exclude 'results/' --exclude 'figures/' --exclude 'checkpoints/' --exclude 'logs/' \
    ./ <user>@<host>:~/dvi/
```

## Bước 2 — Cài môi trường trên server (chỉ lần đầu)

```bash
ssh <user>@<host>
cd ~/dvi
module avail R 2>&1 | head        # xem cluster có module R tên gì
module load R                     # hoặc: conda activate <env có r-base>
Rscript scripts/install_packages.R
```

Nếu R nằm trong conda env chứ không phải `module`, sửa dòng `module load R`
trong `jobs/reproduce/*.slurm` và `jobs/stages/*.slurm` thành
`conda activate <env>`.

## Bước 3 — Submit

Phải `cd` vào thư mục gốc repo rồi mới `sbatch`, vì script dùng
`$SLURM_SUBMIT_DIR` làm `--project-root`.

```bash
cd ~/dvi
sbatch jobs/reproduce/smoke.slurm      # 4 CPU, 16G, 1h  — kiểm tra pipeline
sbatch jobs/reproduce/reduced.slurm    # 8 CPU, 32G, 12h — cấu hình giảm tải
sbatch jobs/reproduce/paper.slurm      # 16 CPU, 64G, 48h — cấu hình bài báo
```

Khi tường thời gian của cluster ngắn hơn thời gian chạy `paper`, tách hai giai
đoạn:

```bash
jid=$(sbatch --parsable jobs/stages/paper-simulations.slurm)
sbatch --dependency=afterok:$jid jobs/stages/paper-analysis.slurm
```

## Bước 4 — Theo dõi

```bash
squeue -u $USER
tail -f logs/reduced-<JOBID>.out
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS
```

Job có `--requeue` và pipeline ghi kết quả tăng dần theo tổ hợp (kịch bản, tương
quan, lần lặp, learner), nên `sbatch` lại đúng script sau khi hết thời gian sẽ
chạy tiếp phần còn thiếu thay vì làm lại từ đầu.

## Bước 5 — Kéo kết quả về máy cá nhân

```bash
rsync -az --partial --info=progress2 \
    <user>@<host>:~/dvi/{results,figures,logs}/ ./
```

---

# Tham khảo

## Tham số của `main.R`

```bash
Rscript main.R --help
```

| Tham số | Mặc định | Ý nghĩa |
| --- | --- | --- |
| `--profile` | `reduced` | `smoke`, `reduced` hoặc `paper` |
| `--threads` | `4` | số luồng BLAS/OpenMP |
| `--stages` | tất cả | danh sách stage, phân tách bằng dấu phẩy |
| `--project-root` | tự dò | thư mục chứa `data/`, `results/`, `figures/` |
| `--energy-data` | `data/ENB2012_data.xlsx` | đường dẫn dữ liệu Energy |
| `--concrete-data` | `data/Concrete_Data.xls` | đường dẫn dữ liệu Concrete |
| `--interval` | `paper` | `paper` (se² = s²/B + c²/n, c = Var(Y)²) hoặc `cross` (se² = s²/B) |
| `--dpi` | `300` | độ phân giải PNG |
| `--language` | `R` | nhãn thư mục đầu ra |

## Chạy từng phần

Các stage: `targets`, `simulations`, `baseline`, `ablation`, `real`, `report`,
`check`. Stage `baseline` phải chạy cùng `simulations` vì bảng so sánh cần kết
quả của mô phỏng chính.

```bash
Rscript main.R --profile paper --stages simulations,baseline
Rscript main.R --profile paper --stages targets,ablation,real,report
```

## Notebook

`src/dvi-r-reproduce.ipynb` dùng lại đúng module trong `R/` và trình bày phần
diễn giải phương pháp. Cần kernel R:

```bash
Rscript scripts/install_packages.R --notebook
Rscript -e 'IRkernel::installspec()'
jupyter notebook src/dvi-r-reproduce.ipynb
```

Mở từ thư mục gốc hay từ `src/` đều được vì notebook tự dò `R/config.R`. Đổi
profile ở cell cấu hình.

## Cấu trúc

```
data/                 ENB2012_data.xlsx (768 x 10), Concrete_Data.xls (1030 x 9)
docs/                 bài báo gốc và yêu cầu đồ án
R/                    module phương pháp
main.R                điểm vào khi chạy bằng Rscript
scripts/              cài đặt môi trường
jobs/reproduce/       script Slurm chạy trọn một profile
jobs/stages/          script Slurm tách giai đoạn cho profile `paper`
src/                  notebook trình bày, dùng lại đúng module trong R/
```

| Module | Nội dung |
| --- | --- |
| `R/config.R` | thư mục gốc, đường dẫn dữ liệu, giới hạn luồng, thứ tự nạp module |
| `R/constants.R` | hằng số, bảng màu, `get_profile`, coverage tham chiếu Table 2 |
| `R/utils.R` | tiện ích số học, chia fold, basis đa thức trực giao, KDE, thống kê an toàn |
| `R/models.R` | nuisance learner: linear, `mgcv::gam`, `grf::regression_forest` |
| `R/estimators.R` | psi_L, psi_0, psi_1, psi_2, psi_3, residual CPI, marginal PFI theo fold |
| `R/inference.R` | tổng hợp fold và khoảng t-Cross |
| `R/simulation.R` | DGP Kịch bản 1–5, checkpoint, mô phỏng chính, benchmark, ablation |
| `R/real_data.R` | nạp, kiểm tra và phân tích Energy và Concrete |
| `R/plots.R` | hình PNG 300 DPI và PDF vector |
| `R/reporting.R` | metadata môi trường, hướng dẫn thực thi, manifest đầu ra |
| `R/pipeline.R` | các stage và hàm điều phối `run_pipeline` |

## Đầu ra

| Đường dẫn | Nội dung |
| --- | --- |
| `results/R/*.csv` | kết quả theo fold, bảng tổng hợp, chẩn đoán, so sánh baseline, ablation, checklist |
| `results/R/environment_<profile>.json` | phiên bản R và package, seed, cấu hình, nền tảng |
| `results/R/output_manifest.csv` | danh sách tệp đầu ra kèm kích thước và thời điểm |
| `figures/R/*.png`, `*.pdf` | hình 300 DPI và bản vector |
| `checkpoints/R/` | dữ liệu mô phỏng, fold assignment, cache tỷ số mật độ |
| `README_EXPERIMENTS.md` | hướng dẫn thực thi sinh tự động kèm SHA-256 dữ liệu |

Các thư mục này được `.gitignore` bỏ qua vì tái tạo được từ mã nguồn.

## Quy ước diễn giải

- Coverage với dưới 30 lần lặp chỉ mang tính mô tả; chỉ profile `paper` (100 lần
  lặp) mới so sánh trực tiếp được với Table 2 của bài báo.
- Khoảng được báo cáo chính là t-Cross theo mục 3.6 của bài báo
  (`se² = s²/B + c²/n`, `c = Var(Y)²`), nên cột `coverage` và
  `coverage_minus_paper` so trực tiếp được với Table 2. Khoảng không hiệu chỉnh
  vẫn được tính và lưu ở các cột `*_cross`; đổi bằng `--interval cross`. Phần
  ablation đánh giá `c = 0`, `c = Var(Y)` và `c = Var(Y)²`.
- Marginal PFI không cùng estimand với tham số decorrelated; vị trí của nó trong
  bảng chỉ mô tả hệ quả của phép nhiễu biên.
