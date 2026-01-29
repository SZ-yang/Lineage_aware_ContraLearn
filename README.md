# Lineage-aware_CL

LCL is a contrastive learning framework that treats inheritable lineage barcodes as a "natural" data augmentation mechanism to isolate subtle fate-determining signals and enable robust representation transfer to unlabeled single-cell datasets.

## Environment requirement

This project requires a Conda environment specified in `environment.yml`.

### Step 1: Download and create the environment

```bash
conda env create -f environment.yml
conda activate LCL
```

## Example usage

### Train the LCL model (semi-supervised)

To see all available arguments:

```bash
python LCL_Main_Semi.py --help
```
Example training command:
```bash
python LCL_Main_Semi.py \
  --inputFilePath <TRAIN_DATA.h5ad> \
  --testFilePath <TEST_DATA.h5ad> \
  --batch_size Batch_size \
  --size_factor 0.04 \
  --unlabeled_per_batch 5 \
  --lambda_penalty 0.05 \
  --temperature 0.5 \
  --output_dir <OUTPUT_DIR> \
  --train_test 1
```
### Feature Extraction
After training, extract embeddings from the base encoder/proj head using the final checkpoint:
To extract embeddings from the base encoder:
```bash
python feature_extraction.py \
  --inputFilePath <INPUT_DATA.h5ad> \
  --batch_size Batch_size \
  --output_dir <OUTPUT_DIR> \
  --resume_from_checkpoint <CHECKPOINT.ckpt> \
  --out_file_name base_embeddings.npy

```
To extract embeddings from the projection head:
```bash
python project_head_extraction.py \
  --inputFilePath <INPUT_DATA.h5ad> \
  --batch_size Batch_size \
  --output_dir <OUTPUT_DIR> \
  --resume_from_checkpoint <CHECKPOINT.ckpt> \
  --out_file_name projection_embeddings.npy
```



## File Organization

```text
main/                                # Python scripts to train and run LCL

analysis/
├── data_processing/                 # Preprocessing and dataset construction
│   ├── cellTag/
│   ├── cellTagMulti/
│   ├── cellTag_cellTagMulti_Integration/
│   └── LARRY/
│
├── LCL/                             # Core LCL experiments and analyses
│   ├── train_test/                  # Train/test split: KNN & KL divergence (Sections 5.2–5.3)
│   ├── train_test_unseen_lineage/   # Supplementary unseen-lineage generalization 
│   ├── simulation/                  # Pseudo-real simulations (Section 5.1)
│   ├── gemli/                       # GEMLI comparisons and memory gene analyses (Section 5.4)
│   ├── gsea/                        # Cospar and GSEA analyses
│   └── plotting/                    # Figure generation
│
├── scVI/                            # scVI baseline experiments
├── supUMAP/                         # Supervised UMAP baseline
├── variancePartition/               # Supplementary variance analyses
│
example_usage.ipynb                  # Minimal example: training LCL on a new dataset
```