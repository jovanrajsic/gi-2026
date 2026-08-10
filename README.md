---

# Single-Cell Analysis of Immune Response to Nanoplastic Particles

**Course:** Genomic Informatics (2026), School of Electrical Engineering, University of Belgrade

## 1. Prerequisites
- **Python:** version 3.12
- **Ubuntu/Debian:** sudo apt install build-essential python3.12-dev

## Environment setup

```bash
# Create and activate environment
python3 -m venv venv
source venv/bin/activate
```

## 3. Data Acquisition
Put files from [Zenodo](https://zenodo.org/records/15866724) in a folder named `data/`.

| Filename | Condition |
| :--- | :--- |
| `filtered_feature_bc_matrix.h5ad` | 40nm PSNPs |
| `filtered_feature_bc_matrix_Sample2.h5ad` | 200nm PSNPs |
| `filtered_feature_bc_matrix_Sample3.h5ad` | 40nm + 200nm Mixture |
| `filtered_feature_bc_matrix_Sample4.h5ad` | Control (No exposure) |

## 4. Run locally

1. Install VSCodium
2. Install Jupyter and Python extensions in VSCodium
3. Set up environment
4. Select kernel from the environment
5. Run

## 5. Outputs
Outputs will be in the `results/` directory.

## Video Presentation
[![Watch the video](https://img.youtube.com/vi/etqQlJw8kg4/0.jpg)](https://www.youtube.com/watch?v=etqQlJw8kg4)

## Disclaimer
I used Google AI Studio for researching code solutions and biological interpretations of the data.
I tried using AI agents like Google Antigravity but it worked out a solution of the project without my interference and I felt that I don't own the solution so I abandoned it.