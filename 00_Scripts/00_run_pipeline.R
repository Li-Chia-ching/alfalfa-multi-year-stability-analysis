# 00_run_pipeline.R
# 主控脚本：一键执行紫花苜蓿项目全流程分析与绘图

cat("=======================================================\n")
cat("Starting Analysis Pipeline...\n")
cat("=======================================================\n")

# 记录开始时间
start_time <- Sys.time()

# 按顺序执行脚本 (假设脚本都在同一目录下)
scripts_to_run <- c(
  "01_import_data.R",
  "02_clean_data.R",
  "03_fit_mixed_model.R",
  "04_estimate_blup.R",
  "05_heritability.R",
  "06_ammi.R",
  "07_gge.R",
  "08_cluster.R",
  "09_export_tables.R",
  "10_make_figures.R"
)

for (script in scripts_to_run) {
  cat(sprintf("\n>>> [%s] Running %s...\n", format(Sys.time(), "%H:%M:%S"), script))
  # 使用 local=new.env() 防止全局环境变量互相污染，保持严格的文件输入输出
  source(script, local = new.env())
}

end_time <- Sys.time()
cat("\n=======================================================\n")
cat("Pipeline completed successfully in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
cat("All Results, CSVs, and Figures have been generated.\n")
cat("=======================================================\n")
