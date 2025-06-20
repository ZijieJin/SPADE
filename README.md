# SPADE

SPADE is a computational pipeline for integration of scRNA-seq data and spatial transcriptomics (ST) data. The main function of SPADE include: 
- Map scRNA-seq data to spatial
- Impute spatial gene expression
- Calculate cell-cell communication scores between spots and detect CCC hotspots 
 
If you have any questions related to scFusion, please visit https://github.com/ZijieJin/SPADE and post them on the *Issues* page or email me: jzj2035198@outlook.com

## Software Prerequisite
SPADE works on R and Python platform with any OS.
- R (tested on 4.2.3)
- R package: stringr, entropy, pracma, RcppML, NMF, matrixStats, progress
- Python (tested on 3.12)
- Python modules: numpy (tested on 1.26.4), pandas (tested on 2.2.3), torch (tested on 2.5.0), scikit-learn (tested on 1.5.2), argparse, os

## System Requirement

To run SPADE properly, Your computer should have:

- 32 GB memory or more
- Nvidia GPU or Apple Silicon GPU

## Data Requirement

SPADE has two mendatory input and other optional inputs:

- Single-cell RNA-seq gene expression data
- Spatial transcriptomic gene expression data
- (optional) Single-cell RNA-seq cell annotations
- (optional) Coordinate of spots
- (optional) Ligand-receptor database
- (optional) downstream regulation database

## Quick Start

Suppose you have prepared the input data `sc.csv` for scRNA-seq expression, `st.csv` for ST gene expression, `location.csv` for ST coordinates, `annotation.csv` for single-cell annotation, run:

`python Run_SPADE.py -s st.csv -r sc.csv -a annotation.csv -c location.csv`

## Full Usage

SPADE have command parameters below:

| Short Argument | Full Argument   | Type    | Default      | Description               |
|:-----|:------------|:---------|:--------------|:--------------------------------------|
|`-s`| `--st`   | `str`   |  Mandatory    | Path to the ST expression file                      |
|`-r`| `--sc`  | `str`   |  Mandatory| Path to the SC expression file                |
|`-a`| `--sc_anno`  | `str`   | ` `          | Path to the SC annotation file                   |
|`-c`| `--st_coord`      | `str` | ` `       | Path to the ST coordinate file             |
|`-o`| `--outputpath` | `str`  | `./`       | Path to save the output file             |
|| `--useAllGenes`    | `boolean`   | `False`          | If set, use all genes instead of marker genes             |

## Commercial Use

For non-academic use, please email Dr. Jin (jzj2035198@outlook.com) to obtain the paid commercial license.

For academic use, source code is licensed under MIT License. 
