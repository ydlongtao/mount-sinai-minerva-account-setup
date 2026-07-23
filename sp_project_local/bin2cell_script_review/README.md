# Cellpose / GEX / bin2cell 脚本审阅包

本目录是从 Minerva 服务器当前项目下载的运行脚本副本，供人工检查。服务器上的脚本和正在运行的任务没有被修改。

## 当前正式修正版主线

1. `scripts/22_cellpose_he_each_sample.py`
   - 读取 2um Visium HD bins 和 H&E 图像。
   - 使用实际源图像 mpp 执行 H&E Cellpose。
   - 生成 `labels_he`，并扩展为 `labels_he_expanded`。
   - 支持通过环境变量测试 `flow_threshold`、`prob_thresh` 和独立输出目录。
2. `scripts/23_bin2cell_raw_counts_each_sample.py`
   - 将 H&E 扩展标签写入 `labels_joint`。
   - 使用 `ov.space.bin2cell()` 创建细胞对象。
   - 从原始 2um H5 重新聚合 raw counts，避免使用已经归一化的矩阵。
   - 输出 cell-level H5AD 和 Space Ranger v4 风格 GeoJSON/H5。
3. `scripts/24_cell_level_qc_marker_each_sample.py`
   - 统计每个细胞的 counts、detected genes 和前列腺癌 marker score。
   - 当前只评分和报告，不自动过滤细胞。
4. `scripts/47_gex_segmentation_compare.py`
   - 用 2um GEX counts 构造表达图像并执行 GEX segmentation。
   - 生成 GEX-only 和 H&E + GEX salvage 两套 cell-level 结果。
   - 绿色区域表示由 GEX 补回的 H&E 未覆盖区域。

## LSF 执行顺序

- `lsf/28_corrected_cellpose_pilot_gpu.lsf`：修正版 Cellpose GPU 试点。
- `lsf/29_corrected_bin2cell_qc_cpu.lsf`：bin2cell 聚合和 QC。
- `lsf/30_corrected_cellpose_cpu.lsf`：无 GPU 时的 CPU Cellpose 重试，使用共享离线模型。
- `lsf/32_cellpose_parameter_sweep_r4_r5.lsf`：R4/R5 的 H&E 参数 sweep。
- `lsf/33_bin2cell_sweep_r4_r5.lsf`：对 6 组 H&E 结果进行 bin2cell 和 overlay。
- `lsf/35_gex_segmentation_r4_r5.lsf`：R4/R5 的 GEX segmentation 和 H&E salvage。
- `lsf/34_bin2cell_sweep_report.lsf`、`lsf/36_gex_segmentation_report.lsf`：生成 HTML 报告。

## 关键数据规则

- 输入表达矩阵使用 `square_002um` 的 raw H5。
- `labels_he_expanded` 是当前 H&E 主标签。
- GEX 结果使用 `labels_gex` 作为次级标签，只填充 H&E 未覆盖 bins。
- `labels_joint` 是最终用于 bin2cell 的联合标签。
- 原始数据目录只读；临时文件和工作文件放在 `/sc/arion/scratch/huangl21/sp_project/`。
- R4/R5 sweep 输出不会覆盖 baseline Cellpose 结果。

## 目录来源

- 服务器项目：`/sc/arion/work/huangl21/sp_project/`
- 本地审阅目录：`sp_project_local/bin2cell_script_review/`
- 下载内容不包含 H5AD、图像或其他大数据文件。
