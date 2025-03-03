> # BGI-Human_plaque
The file "0_QC.R" is a scRNA-seq data quality control script.
The file "1_harmony.R" is a script for scRNA integration, batch effect removal, clustering, and grouping.
The file "2_sterero-seq_annotation.R" is a script that utilizes scRNA annotation data to assist in cell annotation for Stereo-seq.
The file "3_GO&KEGG_enrichment.R" is a script designed for performing GO enrichment and KEGG enrichment analyses using a specific set of genes.
The files "4.1_Trajectory_analysis_monocle2.R", "4.2_Trajectory_analysis_CytoTrace.R", and "4.3_Trajectory_analysis_monocle3.R" are scripts that utilize monocle2, CytoTrace, and monocle3, respectively, to infer the differentiation trajectories of different specific cell subtypes.
The files "5.1_lineage_step1.py", "5.2_lineage_step2.py", and "5.3_lineage_plot.R" are scripts that analyze cell differentiation lineages using mitochondrial SNPs.
The file "6_spatial_bin_size.R" transforms Stereo-seq data into a format with a bin size of 50.
The file "7_random_area.ipynb" is used to randomly select non-overlapping regions with identical contours within the tissue area of a spatial chip.
