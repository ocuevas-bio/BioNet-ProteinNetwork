# 🧬 BioNet

<p align="center">
  <img src="assets/BioNetLogo1.png" alt="BioNet Logo" width="180"/>
</p>

**BioNet** is a desktop application built in Python with a graphical user interface (Tkinter) designed for the structural and network analysis of proteins from PDB files. It enables researchers and students to construct residue contact networks, compute centrality metrics, and interactively visualize 3D molecular structures[cite: 3, 4].

---

## ✨ Key Features
* **3D Network Visualization:** Interactive construction and rendering of amino acid contact networks using Mayavi and NetworkX.
* **Centrality Properties:** Detailed calculation and color-mapped visualization of key network metrics (Degree, Betweenness, Closeness, Clustering Coefficient, Eigenvector, Eigenvector $\lambda_2$, Eccentricity, and Stress).
* **Distance-Based Analysis:** Generation of contact networks using default distance thresholds or custom distance ranges with configurable edge colors.
* **Direct RCSB PDB Fetching:** Load local `.pdb` files or download structures directly via their PDB ID.
* **Data Export:** Export calculated centrality metrics, node classifications, and network summaries to Excel (`.xlsx`), CSV, JSON, or TXT formats[cite: 3, 4].

---

## 🛠️ Technologies & Libraries
* **Python 3**
* **Tkinter** (Graphical User Interface)
* **NetworkX** (Graph theory and network metrics)
* **BioPython** (PDB structure parsing)
* **Mayavi & Matplotlib** (3D interactive rendering and 2D plotting)
* **Pandas & NumPy** (Data manipulation and numerical computing)[cite: 3, 4, 6]

---

## 📂 Project Structure
```text
BioNet/
├── assets/
│   ├── BIONET.ico
│   └── BioNetLogo1.png
├── notebooks/
│   └── BioNet1_7.ipynb
├── BioNet.py
├── requirements.txt
└── README.md
```[cite: 1, 2, 3, 4, 5]

---

## 🚀 Installation & Usage

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/your-username/BioNet.git](https://github.com/your-username/BioNet.git)
   cd BioNet
## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

[![DOI](https://zenodo.org/badge/1456117046.svg)](https://doi.org/10.5281/zenodo.24203599)
