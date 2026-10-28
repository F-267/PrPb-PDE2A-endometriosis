# Analysis code for the manuscript

Code for the computational analyses in Figure 1, Fig. 2I and Fig. S1:
network toxicology + machine learning (LASSO/RF/SVM/Venn), nomogram,
multi-classifier ROC, SHAP interpretation, KEGG GSEA, and mouse GSEA.

## Scripts

| Script | Purpose | Figure | Dataset |
| --- | --- | --- | --- |
| `01_LASSO_RF_SVM_Venn_with_CV.R` | LASSO, random forest, SVM permutation importance, three-way Venn intersection; 10-fold CV accuracy/error curves | Fig. 1D-H | GSE141549 |
| `02_nomogram.R` | Nomogram, ROC, calibration and DCA curves | Fig. 1I | GSE51981 |
| `03_ROC_and_SHAP.R` | Ten-classifier ROC comparison and SHAP interpretation | Fig. 1J, Fig. S1 | GSE141549 |
| `04_GSEA_KEGG.R` | KEGG GSEA of high- vs low-PDE2A groups | Fig. 1L | GSE25628 |
| `05_GSEA_mouse_Pde2a.R` | GSEA of high- vs low-Pde2a mouse lesions (MSigDB C2:CP) | Fig. 2I | mouse RNA-seq (EM_vs_NC) |

## Required input files

Each script reads its inputs from the working directory set by `setwd(...)`:

- `01_LASSO_RF_SVM_Venn_with_CV.R`: `02normal_exp_uniq_gene.csv` (GSE141549), `cogene.csv`, `group.csv`
- `02_nomogram.R`: `02normal_exp_uniq_gene.csv` (GSE51981), `group.csv`, `rf_lasso_svm_gene.csv`
- `03_ROC_and_SHAP.R`: `merge.normalize.txt` (GSE141549), `interGenes.txt` (14 hub genes), `refer.methodLists.txt` (classifier list)
- `04_GSEA_KEGG.R`: `merge.normalize.txt` (GSE25628), `c2.cp.kegg.Hs.symbols.gmt`
- `05_GSEA_mouse_Pde2a.R`: `EM_vs_NC.normalize.txt` (mouse), `m2.cp.v2026.1.Mm.symbols.gmt`

## Data

The expression matrices are public GEO datasets:
- GSE141549 (machine-learning / ROC)
- GSE51981 (nomogram)
- GSE25628 (GSEA; also used for the DEG step with GSE7305)
- GSE7305 (DEG)

They are not bundled because of size. Small gene lists are included in `data/`:
- `cogene.csv` - candidate genes (PrPb targets intersected with EM DEGs)
- `rf_lasso_svm_gene.csv` - the 14 hub genes used for nomogram and ROC
- `interGenes_14.txt` - the same 14 hub genes as a plain list
- `EM_vs_NC.normalize.txt` - mouse lesion expression matrix (6 samples: 3 control, 3 EM) for the mouse GSEA (Fig. 2I)
- `m2.cp.v2026.1.Mm.symbols.gmt` - mouse MSigDB C2:CP gene sets for the mouse GSEA

## Notes

- R version 4.5.x. Install missing packages with `install.packages()` or
  `BiocManager::install()`.
- Each script contains a `setwd(...)` line; change it to your local directory
  before running.
- The 14 hub genes were derived from GSE141549 and used for the nomogram
  (GSE51981), ROC (GSE141549) and GSEA (GSE25628).
