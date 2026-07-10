# Data Structure Notes

Raw data root:

```bash
/sc/arion/projects/DiseaseGeneCell/Huang_lab_data/SpatialTranscriptome
```

Observed 10x Space Ranger output pattern:

```text
<group>/<sample>/outs/
  binned_outputs/
    square_002um/
      filtered_feature_bc_matrix.h5
    square_008um/
      filtered_feature_bc_matrix.h5
      spatial/
    square_016um/
      filtered_feature_bc_matrix.h5
      spatial/
  segmented_outputs/
    filtered_feature_cell_matrix.h5
  spatial/
  metrics_summary.csv
  web_summary.html
  feature_slice.h5
  barcode_mappings.parquet
  possorted_genome_bam.bam
  molecule_info.h5
  *.cloupe
```

Version 1 reads:

- `binned_outputs/square_008um/filtered_feature_bc_matrix.h5`
- `binned_outputs/square_008um/spatial/`
- `binned_outputs/square_016um/filtered_feature_bc_matrix.h5` for manifest and QC shape only
- `segmented_outputs/filtered_feature_cell_matrix.h5`
- `metrics_summary.csv`

Version 1 does not read:

- `square_002um`
- `possorted_genome_bam.bam`
- `molecule_info.h5`
- `*.cloupe`

The local structure snapshot retained outside this package is:

```bash
/Users/huangfulongtao/Documents/Minerva服务器账户配置/SpatialTranscriptome_structure_20260709_110442.md
```
