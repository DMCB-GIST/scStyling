# 🚀 Running Commands in each model Folder

This repository provides several scripts to run different models.  
All commands should be executed inside the each model folder.  

---

## 📌 Available Scripts

| Model                          | Script Location                                                                 |
|--------------------------------|---------------------------------------------------------------------------------|
| **scStyling**                  | `run_scStyling.sh`                                                              |
| **Simple Classifier**          | `run_Simple_classifier.sh`                                                      |
| **SCAD**                       | `SCAD/model/run_SCAD.sh`                                                        |
| **SCAD (with foundation embeddings)** | `SCAD/model/run_SCAD.sh` <br> `scFoundation/SCAD/model/run_SCAD_with_embeddings.sh` |
| **scDEAL (bulkmodel)**         | `scDEAL/run_bulkmodel.sh`                                                       |
| **scDEAL (scmodel)**           | `scDEAL/run_scmodel.sh`                                                         |

---

## ⚙️ Command Arguments

All scripts take the following arguments:

- **`data_type`**  
  - Options: `Bulk`, `scStyleBulk`  
  - Specifies the type of input data.

- **`foundation_model`**  
  - Options: `gexpr`, `scRobust`, `scGPT`, `scFoundation`  
  - Defines whether to use raw gene expression (`gexpr`) or foundation model embeddings.

---

## Download Data
https://zenodo.org/records/17287489


## 🗂️ Dataset Tree (Markdown Links)

- [data/](./data)  
  - [embeddings/](./data/embeddings)  
  - [gexpr/](./data/gexpr)  

- [SCAD/](./SCAD)  
  - [data/](./SCAD/data)  

- [scFoundation/](./scFoundation)  
  - [SCAD/](./scFoundation/SCAD)  
    - [data/](./scFoundation/SCAD/data)

---

## ⚙️ Usage Example

```bash
# Clone the repository
git clone https://github.com/your-username/scStyling.git
cd scStyling

# Explore datasets
ls data/
ls SCAD/data/
ls scFoundation/SCAD/data/
