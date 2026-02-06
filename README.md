# Nanovaccines
README.md
LICENSE
.gitignore
renv.lock (或 sessionInfo.txt)

code/                 # 入口脚本（最重要）
  00_setup.R
  02_make_figures.R
  03_make_supp_figs.R
  10_paths.R          # 统一路径（强烈推荐）

R/                    # 共享函数（从各子项目抽出来）
  data_prep.R
  data_combine.R
  plot_combined.R
  plot_separate.R
  plot_pca.R
  ...

data/                 # “最终用于画图”的数据（尽量统一）
  raw/                # 原始数据（若可公开）
  processed/          # 预处理后（csv/rds/rda）
  README_data.md

outputs/
  figures/
  tables/
  logs/

archive/              # 你 zip 里原来的每个文件夹原样放这
 OCR/
  code_invivo/
  ...
