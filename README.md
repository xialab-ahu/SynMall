# 1 🌐 SynMall
> **SynMall** is a comprehensive database dedicated to human whole-genome synonymous mutations, integrating large-scale annotation, feature extraction, and literature-based evidence.  
> The web portal is available here: [**SynMall**](https://bioinfo.ahu.edu.cn/synMall/#/home)


<img width="2498" height="1223" alt="image" src="https://github.com/user-attachments/assets/e84b7ed0-e4ad-429d-b157-de8bbd62dd6a" />


This repository provides all the code, figures, and results used in our manuscript to ensure full reproducibility.  
To reproduce the work, please follow the steps below.


# 2 🚀 Environment Setup

```bash
git clone https://github.com/ToolForVol/SynMall.git
cd SynMall
conda env create -f environment.yml
conda activate SynMall
```

# 3 📦 Data Download

Run the script below to download all required files for analysis:

```bash
cd SynMall
bash ./download.sh
```

# 4 ⚙️ Reproduction Workflow

1. `./Code/1-SynMall-Prediction.ipynb` - Construction, training, prediction, and benchmarking of the synScore model.
2. `./Code/2-Whole-Genome-Prediction.py` - Genome-wide prediction using synScore.
3. `./Code/3-SynScore-Analysis.ipynb` - Statistical analysis and visualization based on synScore.
4. `./Code/4-Feature-Module-Visualization.ipynb` - Visualization code in our manuscript `Figure 2`.
5. `Pre-computed genome-wide (hg38) synScore predictions` can also be obtained via `./download.sh`.

# 5 🧬 Citation

```cite
In preparation.
```
