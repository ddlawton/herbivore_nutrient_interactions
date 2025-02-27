# Exploring Nutrient Availability and Herbivorous Insect Population Dynamics Across Multiple Scales

## Overview
This repository contains the data and code supporting the publication:

**Lawton et al. (2024).** *Exploring Nutrient Availability and Herbivorous Insect Population Dynamics Across Multiple Scales.* Oikos. DOI: _[To be added]_

The study examines how nutrient availability influences herbivore populations at different spatial scales, using a combination of field experiments, remote sensing data, and spatial modeling.

## Repository Structure
The repository is organized into major sections, each containing analysis scripts and associated data:

- **`field_populations/`** - Analyses related to field population studies.
  - `choice_intake_target_experiment/`
  - `no_choice_intake_target_experiment/`
- **`field_cages/`** - Investigations of locust survival, growth, and nutrient intake.
  - `grass_nutrient_analysis/`
  - `locust_survival_growth_analysis/`
  - `locust_intake_target_analysis/`
- **`spatial_modeling/`** - Spatial analyses and predictive modeling.
  - `fishnet_grid_construction/`
  - `remote_sensing_data_extraction/`
  - `raw_data_visualization/`
  - `modeling/`
  - `post_modeling_visualization/`

Notebooks within each subfolder are numbered sequentially (e.g., `01_analysis.Rmd`, `02_modeling.Rmd`) to indicate the order of execution.

## Requirements
This project uses both R and Python for data analysis. Dependencies are managed via `renv` (R) and a virtual environment (`myenv`) for Python.

### Software Requirements
- **Python 3.x** (dependencies managed via `myenv`)
- **R 4.x** (dependencies managed via `renv`)

### Installation
Clone the repository and set up the required environments:

```bash
git clone https://github.com/username/repository-name.git
cd repository-name
```

#### Python Setup
Activate the Python virtual environment:
```bash
conda activate myenv  # or use `source myenv/bin/activate` if using virtualenv
```

#### R Setup
Restore the R package environment:
```r
install.packages("renv")
renv::restore()
```

## Data Availability
The datasets used in this study are archived on Zenodo and will be available upon publication. DOI: _[To be added]_

## Reproducibility
To reproduce the analyses, follow these steps:
1. Ensure all dependencies are installed (see above).
2. Navigate to the relevant analysis subfolder (e.g., `spatial_modeling/modeling/`).
3. Run the scripts in sequential order.

## Citation
If you use this dataset or code, please cite:

**Lawton et al. (2024).** *Exploring Nutrient Availability and Herbivorous Insect Population Dynamics Across Multiple Scales.* Oikos. DOI: _[To be added]_

## License
This repository is licensed under the MIT License. See `LICENSE` for details.

## Contact
For questions or feedback, please contact: _ddlawton@asu.edu_.

