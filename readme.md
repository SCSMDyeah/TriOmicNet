# TriOmicNet

TriOmicNet is a multi-layer network diffusion framework for cancer driver-gene identification. It integrates three complementary components derived from gene expression, somatic mutation, miRNA–gene regulation, and multiple protein interaction networks:

1. control ability score;
2. regulatory potential score; and
3. multi-network diffusion score.

The three component scores are standardized and integrated by singular value decomposition (SVD) to produce the final `TriOmicNetScore`.

## Workflow

The executable workflow contains five steps:

1. expression data preprocessing;
2. multilayer network construction;
3. calculation of the control ability score and regulatory potential score;
4. calculation of the multi-network diffusion score; and
5. SVD-based integration of the three component scores.

All five steps are executed by `code/TriOmicNet.R`. The multilayer-network construction step is not a descriptive placeholder: the main script explicitly loads `code/construct_layer.R` and calls `construct_layer()`.

## Repository structure

```text
TriOmicNet/
├── code/
│   ├── TriOmicNet.R
│   ├── DE_Score.R
│   ├── construct_layer.R
│   ├── random_walk.R
│   └── remrna_exp.m
├── data/
│   ├── BRCA/
│   │   ├── brca_gene_exp.RData
│   │   ├── brca_mut.RData
│   │   ├── brca_cancer_si.rda
│   │   ├── brca_normal_si.rda
│   │   └── remrna_exp_brca.txt        # processed file; may also be placed in code/
│   ├── LUAD/
│   │   ├── luad_gene_exp.RData
│   │   ├── luad_mut.RData
│   │   ├── luad_cancer_si.rda
│   │   ├── luad_normal_si.rda
│   │   └── remrna_exp_luad.txt        # required to run LUAD
│   ├── PRAD/
│   │   ├── prad_gene_exp.RData
│   │   ├── prad_mut.RData
│   │   ├── prad_cancer_si.rda
│   │   ├── prad_normal_si.rda
│   │   └── remrna_exp_prad.txt        # required to run PRAD
│   ├── network.rda
│   ├── protein_info.rda
│   ├── mir_gene_network.rda
│   ├── HINT.RData
│   ├── CPDB.RData
│   └── MULTINET.RData
└── results/
    ├── BRCA/
    ├── LUAD/
    └── PRAD/
```

The main script also accepts processed expression files located in `code/`, `data/`, or the repository root. Both `.txt` and `.txt.gz` formats are supported.

## Software environment

The workflow was developed in R. The required packages are:

```r
install.packages(c("readr", "igraph", "irlba"))
```

The script checks these packages before starting and reports a clear error if one is missing.

## Input data

For each cancer type, the following four cancer-specific R data files are required:

```text
data/<CANCER>/<prefix>_gene_exp.RData
data/<CANCER>/<prefix>_mut.RData
data/<CANCER>/<prefix>_cancer_si.rda
data/<CANCER>/<prefix>_normal_si.rda
```

where `<CANCER>` is `BRCA`, `LUAD`, or `PRAD`, and `<prefix>` is the corresponding lowercase name.

The processed expression matrix is expected to be named:

```text
remrna_exp_brca.txt
remrna_exp_luad.txt
remrna_exp_prad.txt
```

`code/remrna_exp.m` can be used to generate these files. A processed BRCA file may be supplied so that the BRCA example can be run without MATLAB. LUAD and PRAD require their corresponding processed expression files before those workflows can be executed.

## How to run

Clone or download the repository and open a terminal in the repository root.

### Run the default BRCA example

```bash
Rscript code/TriOmicNet.R
```

This is equivalent to:

```bash
Rscript code/TriOmicNet.R BRCA
```

### Run LUAD

```bash
Rscript code/TriOmicNet.R LUAD
```

### Run PRAD

```bash
Rscript code/TriOmicNet.R PRAD
```

### Run all three cancer types

```bash
Rscript code/TriOmicNet.R ALL
```

The `ALL` command runs BRCA, LUAD, and PRAD sequentially. All three processed expression files must be present.

## Explicit execution of the five analysis steps

### Step 1: Expression data preprocessing

The main script loads the cancer-specific expression data, calculates the mean, variance, standard deviation, fluctuation function, and expression threshold, and imports the processed expression matrix produced by `remrna_exp.m`.

Normal expression samples are identified from sample names containing `11A` or `11B`. The remaining samples are treated as cancer samples. Differential expression scores are then calculated with:

```r
source(file.path(code_dir, "DE_Score.R"))
DE_05 <- DE_Score(normal_s, cancer_s, 0.5)
mirna_value <- DE_Score(normal_si, cancer_si, 0.5)
```

### Step 2: Multilayer network construction

No separate manual command is required. The following commands are executed automatically inside `code/TriOmicNet.R`:

```r
source(file.path(code_dir, "construct_layer.R"))
construct_layer(mut_data, network, protein_info, DE_05)
```

This step uses the mutation matrix, shared biological network, protein information, and differential expression scores. The main script then locates the generated files required by the following steps.

The standard generated file names are:

```text
layer_2_<prefix>.csv
C1_<prefix>.csv
C2_<prefix>.csv
STRING_adj_<prefix>.csv
```

For compatibility with the original implementation, the main script also recognizes the legacy names:

```text
第二层所有节点集合_<prefix>.csv
STRING要构成邻接矩阵的形式_<prefix>.csv
```

The generated files may be written to the repository root, `data/`, `data/<CANCER>/`, `code/`, or `results/<CANCER>/`; the main script checks these locations automatically.

### Step 3: Control ability and regulatory potential scores

The control ability score is calculated by combining the differential miRNA scores with the miRNA–gene regulatory network.

The regulatory potential score is calculated from the direct and indirect regulatory scores:

```r
regulatory_potential <- pmax(S_Dir, S_Ind)
```

### Step 4: Multi-network diffusion score

Random walk with restart is performed on four biological interaction networks:

- STRING;
- HINT;
- CPDB; and
- MULTINET.

The revised script loads CPDB and MULTINET in separate environments. This is necessary because the original `.RData` files may both store their adjacency matrix under the internal object name `PPI`; separate loading prevents one network from overwriting another.

The four network scores are integrated using their mean score and a cross-network consistency coefficient:

```r
consistency_score <- min_score / max_score
multi_network_score <- mean_score * consistency_score
```

The TNI calculation uses all STRING edges and the actual network size. It no longer uses a hard-coded value such as `7686`, and the corrected implementation writes the result to `TNI` rather than to an undefined `ECC` object.

### Step 5: SVD integration

The control ability score, regulatory potential score, and multi-network diffusion score are standardized by Z-score normalization and integrated by SVD:

```r
S <- cbind(
  zscore(control_ability_score),
  zscore(multi_network_diffusion_score),
  zscore(regulatory_potential_score)
)

res <- irlba::irlba(S, nv = 1, nu = 1)
TriOmicNetScore <- abs(res$u[, 1] * res$d[1])
```

If `irlba()` cannot be used for a particular input matrix, the script automatically falls back to the base R `svd()` implementation.

## Output

Results are written to a cancer-specific directory:

```text
results/BRCA/final_score.csv
results/LUAD/final_score.csv
results/PRAD/final_score.csv
```

`final_score.csv` is sorted in descending order of `TriOmicNetScore` and contains a `Rank` column, the final score, and all component scores used by the workflow.

The expression-preprocessing intermediate files are also written to the same cancer-specific result directory:

```text
gene_exp1_<prefix>.txt
mrna_exp_<prefix>.csv
means_<prefix>.txt
```

## Reproducing the main ranking results

To reproduce the three cancer-specific driver-gene rankings reported in the manuscript, run:

```bash
Rscript code/TriOmicNet.R ALL
```

Alternatively, run the three analyses separately:

```bash
Rscript code/TriOmicNet.R BRCA
Rscript code/TriOmicNet.R LUAD
Rscript code/TriOmicNet.R PRAD
```

Each command executes the complete five-step core workflow, including multilayer network construction through `construct_layer()`.


## Notes

The script determines the repository root from the location of `code/TriOmicNet.R`, so it is not dependent on absolute Windows paths. Nevertheless, running the commands from the repository root is recommended.
