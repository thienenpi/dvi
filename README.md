# Tái thực nghiệm Decorrelated Variable Importance

Tái thực nghiệm phương pháp trong Verdinelli và Wasserman, *Decorrelated Variable
Importance* (JMLR 25, 2024 — bản PDF trong `docs/`): năm kịch bản mô phỏng của
bài báo, ba nuisance learner, thủ tục suy luận t-Cross, hai baseline đối chứng
(residual CPI, marginal PFI) và phân tích hai bộ dữ liệu thực Energy Efficiency
và Concrete Compressive Strength.

## Cấu trúc

```
data/                 ENB2012_data.xlsx, Concrete_Data.xls
docs/                 bài báo gốc
R/                    module phương pháp
main.R                điểm vào khi chạy bằng Rscript
jobs/reproduce/       script Slurm chạy trọn một profile
jobs/stages/          script Slurm tách giai đoạn cho profile `paper`
src/                  notebook trình bày, dùng lại đúng module trong R/
results/, figures/, checkpoints/, logs/    đầu ra (sinh khi chạy)
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

## Cài đặt môi trường

Yêu cầu R từ 4.2 trở lên. Repository đã kèm sẵn hai tệp dữ liệu trong `data/`
nên không cần tải thêm.

Cài R (chọn một cách):

```bash
sudo apt-get install -y r-base            # Debian/Ubuntu
conda create -n dvi -c conda-forge r-base # conda
brew install r                            # macOS
```

Cài package R (`jsonlite`, `readxl`, `mgcv`, `grf`; thêm `--notebook` để cài
`IRdisplay` và `IRkernel` cho notebook):

```bash
Rscript scripts/install_packages.R
Rscript scripts/install_packages.R --notebook
```

Kiểm tra nhanh toàn bộ pipeline (vài phút):

```bash
Rscript main.R --profile smoke --threads 4
```

## Chạy thực nghiệm

```bash
Rscript main.R --profile smoke   --threads 4     # kiểm tra nhanh toàn bộ đường đi
Rscript main.R --profile reduced --threads 8     # cấu hình giảm tải
Rscript main.R --profile paper   --threads 16    # cấu hình của bài báo
Rscript main.R --help
```

Tham số: `--profile`, `--threads`, `--dpi`, `--project-root`, `--energy-data`,
`--concrete-data`, `--stages`, `--language`.

Chạy từng phần bằng `--stages` với danh sách trong
`targets,simulations,baseline,ablation,real,report,check`; `baseline` phải chạy
cùng `simulations`.

```bash
Rscript main.R --profile paper --stages simulations,baseline
Rscript main.R --profile paper --stages targets,ablation,real,report
```

## Notebook

`src/dvi-r-reproduce.ipynb` dùng lại đúng module trong `R/` và trình bày phần
diễn giải phương pháp. Cần kernel R (`IRkernel`); mở từ thư mục gốc hoặc từ
`src/` đều được vì notebook tự dò `R/config.R`. Đổi profile ở cell cấu hình.

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

## Slurm

```bash
sbatch jobs/reproduce/smoke.slurm
sbatch jobs/reproduce/reduced.slurm
sbatch jobs/reproduce/paper.slurm
```

Job phải được nộp từ thư mục gốc của repository vì script dùng
`$SLURM_SUBMIT_DIR` làm `--project-root`. Log ghi vào `logs/`.

Mô phỏng chính lưu kết quả tăng dần theo tổ hợp (kịch bản, tương quan, lần lặp,
learner) trong `results/R/simulation_<profile>.csv`, nên nộp lại đúng script sau
khi hết thời gian sẽ chạy tiếp phần còn thiếu. Với cluster có tường thời gian
ngắn, dùng `jobs/stages/paper-simulations.slurm` rồi
`jobs/stages/paper-analysis.slurm`.

## Dữ liệu

`data/ENB2012_data.xlsx` (768 × 10) và `data/Concrete_Data.xls` (1030 × 9).
Pipeline kiểm tra kích thước và ghi SHA-256 vào metadata trước khi chạy. Khi dữ
liệu nằm nơi khác, dùng `--energy-data` và `--concrete-data`, hoặc biến môi
trường `DVI_ENERGY_DATA_PATH` và `DVI_CONCRETE_DATA_PATH`.

## Quy ước diễn giải

- Coverage với dưới 30 lần lặp chỉ mang tính mô tả; chỉ profile `paper` (100 lần
  lặp) mới so sánh trực tiếp được với Table 2 của bài báo.
- Khoảng được báo cáo chính là t-Cross không có số hạng hiệu chỉnh
  (`se² = s²/B`). Khoảng theo `se² = s²/B + c²/n` vẫn được tính và lưu; phần
  ablation đánh giá `c = 0`, `c = Var(Y)` và `c = Var(Y)²`.
- Marginal PFI không cùng estimand với tham số decorrelated; vị trí của nó trong
  bảng chỉ mô tả hệ quả của phép nhiễu biên.
