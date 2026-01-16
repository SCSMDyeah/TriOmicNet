#**TriOmicNet: A Novel Multi-layer Network Diffusion Approach Integrating Three Types of Scores to Identify Cancer Driver Genes**

You can directly run **TriOmicNet.R**, as the data import from the `data` folder, along with the reading of `DE_Score.R`, `construct_layer.R`, and `random_walk.R`, has already been written in the code. The following 5 steps are the main steps of **TriOmicNet**, and you can refer to **Figure 1** in the manuscript for better understanding:

**Step 1:** Data processing. The file for filtering the expression matrix is **remrna_exp.m**. If you don't have MATLAB at the moment, you can directly use **remrna_exp_brca.txt.gz** (this file has already been processed for BRCA cancer).

**Step 2:** Construct the multilayer network.

**Step 3:** Calculate control ability score and regulatory potential score.

**Step 4:** Calculate multi-network score.

**Step 5:** Calculate the final score.
