Below is a sample **README.md** that guides users through your QTL analysis workflow and clarifies the purpose of each script step. It also suggests places to embed images for improved clarity. Feel free to adapt it to your specific context, folder structure, or naming conventions.

---

# README

## Overview

This repository contains R scripts and supporting files for **1-LOD Interval QTL Mapping** using the [R/qtl](https://rqtl.org/) package. The workflow includes:
1. **Data Import and Cleaning**
2. **Genotype Probability Calculation**
3. **Initial QTL Scans with Permutations**
4. **QTL Model Refinement and 1-LOD Interval Extraction**
5. **Effect Plotting and Summaries**

The primary script is [`qtl_analysis.R`](./qtl_analysis.R) (or whatever your file is called), which walks through the entire process from reading your cross data to generating final QTL models and diagnostic plots.

---

## Requirements

- **R** (version 4.0 or higher recommended)
- **R packages:**
  - `qtl`  
  - `rlecuyer`  
  - `snow`  
- A CSV file containing your genotypic and phenotypic data (e.g., `Mapping181_11.csv`)

> **Tip**: Make sure your CSV is formatted as expected by R/qtl (markers in columns, individuals in rows, etc.). You can find more details in the [R/qtl sample data documentation](https://rqtl.org/docsforbroman/qtl/html/sampledata.html).

---

## Project Structure

A typical directory layout might look like:

```
.
├── qtl_analysis.R             # Main R script
├── Mapping181_11.csv          # Example CSV with genotype/phenotype data
├── README.md                  # This file
├── /images                    # Folder for embedding images into the README
└── <other files...>
```

---

## Usage

1. **Open R or RStudio**  
2. **Install required packages** (if you haven’t already):  
   ```r
   install.packages("qtl")
   install.packages("rlecuyer")
   install.packages("snow")
   ```
3. **Set your working directory** to this project folder, or use `setwd()` as indicated in the script.
4. **Run the script** in your R environment:
   ```r
   source("qtl_analysis.R")
   ```
5. The script will:
   - Read in your CSV data.
   - Optionally calculate the best genotyping error probability (commented out by default).
   - Check for normality and batch effects (commented out by default).
   - Perform the genome scans (`scanone`) and permutations.
   - Generate PDF plots of QTL profiles and genotype images.
   - Export text files summarizing QTL hits, intervals, LOD scores, etc.

**All outputs** (PDFs, text summaries, and CSV files) will be saved in the same working directory, unless you customize the file paths.

---

## Step-by-Step Outline

Below is a concise outline of what happens in the script, along with notes on which outputs are generated at each step:

1. **Housekeeping**  
   - Sets working directory, loads packages.

2. **Data Import**  
   - Reads CSV with genotypic/phenotypic data using `read.cross()`.
   - Applies `jittermap()` to slightly shift overlapping marker positions.

3. **Genotype Probabilities**  
   - Sets `error.prob = 0.025` and runs `calc.genoprob()` to compute genotype probabilities.

4. **(Optional) Normality & Batch Effect Diagnostics**  
   - Checks normality with Shapiro-Wilk (commented out by default).  
   - Plots batch effects in `Batch effects.pdf` if you uncomment that section.

5. **Graphical Genotype & Genetic Map**  
   - Saves genotype heatmap to `Graphical Genotype.pdf`.  
   - Saves genetic map overview to `Estimated Genetic Map.pdf`.

6. **Initial QTL Scan**  
   - Uses Haley-Knott regression (`method="hk"`) in `scanone()`.  
   - Performs 1000 permutations to derive significance thresholds.  
   - Summaries are saved in:
     - `10_Initial QTL hits by phenotype.txt`  
     - `10_Summary of top hits by phenotype.txt`  
     - `11_LOD scores for every marker.txt`  

7. **QTL Profiles & Threshold Lines**  
   - Plots genome-wide LOD curves in `QTL Plots.pdf` with threshold lines at \(\alpha = 0.05\).

8. **Marker Map Positions**  
   - Writes `13_Genetic Map Positions.csv` with all marker names and positions.

9. **Genotype Simulation**  
   - Calls `sim.geno()` for 500 genotype draws, supporting effect plotting later.

10. **Model Refinement** (Loop Over Phenotypes)  
    - **Add QTL**: `addqtl()` searches for additional peaks conditional on existing QTL.  
    - **Stepwise Selection**: `stepwiseqtl()` builds the best additive QTL model.  
    - **Refine QTL**: `refineqtl()` polishes positions.  
    - **Export**:
      - `14_Additional QTL hits by phenotype.txt`
      - `15_Additional QTL hits from stepwise analysis.txt`
      - `16_Summary of Final QTL Intervals.txt`  
      - `17_ANOVA results and QTL effect estimates.txt`  
      - **Effect Plots**: One PDF per QTL, e.g. `18_Marker effect plots * .pdf`
      - **Means & SE**:  `19_means and SE.txt`

11. **No QTL Case**  
    - If a phenotype has no peaks above threshold, the script appends a “no QTL found” message in the relevant text files.

---

## Suggested Places to Insert Images

To make this README more visually informative, you might consider embedding sample plots or flowcharts. For instance:

1. **Flowchart of the Analysis Pipeline**  
   - Illustrate the sequence from `read.cross()` → `calc.genoprob()` → `scanone()` → `addqtl()` → `stepwiseqtl()` → `refineqtl()`.
   - Place this near the “Step-by-Step Outline” section.

   ```markdown
   ![Workflow Flowchart](images/workflow_flowchart.png)
   ```

2. **Example QTL Profile**  
   - After you generate `QTL Plots.pdf`, consider saving one representative plot as a PNG and embedding it in your README to show what the LOD peaks look like.

   ```markdown
   ![Example LOD Profile](images/example_lod_profile.png)
   ```

3. **Genotype Heatmap**  
   - Similarly, you can take a screenshot or save a PNG from the `Graphical Genotype.pdf` to illustrate the genotype patterns across individuals and markers.

   ```markdown
   ![Genotype Heatmap](images/genotype_heatmap.png)
   ```

To embed images in Markdown, place the relevant `.png/.jpg/.gif` inside an `images` folder (or any folder you like), and reference it like so:

```markdown
![Descriptive Alt Text](images/your_image_filename.png)
```

---

## Interpreting Key Outputs

1. **`10_Initial QTL hits by phenotype.txt`**  
   - Lists the top markers (and approximate positions) exceeding the significance threshold in the initial scan.

2. **`QTL Plots.pdf`**  
   - Contains the LOD curve for each phenotype. Use these curves to visually inspect QTL peaks and compare LOD thresholds.

3. **`16_Summary of Final QTL Intervals.txt`**  
   - Reports the flanking markers that define a 1-LOD interval around each QTL peak.

4. **`17_ANOVA results and QTL effect estimates.txt`**  
   - Provides model-based inference (ANOVA) for QTL effects and the proportion of variance explained.

5. **Effect Plots** (file names typically `18_Marker effect plots <PHENO> Q<index>.pdf`)  
   - Show how each genotype class differs for the trait, clarifying the effect of a given QTL.

---

## Troubleshooting

- **File Permissions**: Ensure you have write access to the output directory.
- **Missing Markers or Phenotypes**: Double-check your CSV columns match the format R/qtl expects.
- **No QTL Found**: It might be normal if the trait has complex inheritance or not enough statistical power.

---

## License

Include your preferred license here (e.g., MIT, GPL-3.0, etc.). If you have no preference, you can remove this section or simply say “All rights reserved.”

---

## Contact

For questions or issues with this analysis, please contact:
- [Your Name & Email]

Or open an issue in this Git repository (if hosted on GitHub, GitLab, etc.).

---

*End of README*
