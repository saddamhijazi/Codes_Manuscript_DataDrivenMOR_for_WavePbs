# Code for Data-Driven Model Order Reduction for Wave Propagation in Materials with Damage and Nonlinearities

This repository contains the numerical source code used to generate the results presented in the manuscript:

> **Data-driven model order reduction for wave propagation in materials with damage and nonlinearities**

by Saddam Hijazi, Nikiema Fulgence, Hannah Burmester, Natalie Rauter, and Carmen Gräßle.

---

## Associated Dataset

The dataset required to reproduce the numerical results presented in the manuscript is archived on Zenodo:

- **Zenodo DOI:** [10.5281/zenodo.23194087](https://doi.org/10.5281/zenodo.23194087)

---

## Code Structure

The repository is organized according to the numerical examples presented in the manuscript:

```text
Codes_Manuscript_DataDrivenMOR_for_WavePbs/
├── NumericalExample4.1/
│   ├── ...
│   └── ...
├── NumericalExample4.2/
│   ├── ...
│   └── ...
├── NumericalExample4.3/
│   ├── ...
│   └── ...
└── Matlab_Routines/
```

---

## Software Requirements

The code was executed and tested using the following software environment:

- **MATLAB:** R2023b
- **Python:** 3.10.12
- **TensorFlow:** 2.18.0
- **NumPy:** 1.26.4
- **SciPy:** 1.12.0

For complete details on the software environment, installation, and step-by-step reproduction procedure, please refer to [`ReproductionInstructions.pdf`](./ReproductionInstructions.pdf).

---

## Reproducibility Workflow

To reproduce the numerical experiments:

1. **Clone this repository:**

   ```bash
   git clone https://github.com/saddamhijazi/Codes_Manuscript_DataDrivenMOR_for_WavePbs.git
   cd Codes_Manuscript_DataDrivenMOR_for_WavePbs
   ```

2. **Download the dataset** from [Zenodo](https://doi.org/10.5281/zenodo.23194087).

3. **Extract the dataset** into the directory expected by the corresponding numerical example and reproduction instructions.

4. **Consult [`ReproductionInstructions.pdf`](./ReproductionInstructions.pdf)** for the required software environment, execution parameters, and step-by-step instructions.

5. **Run the scripts** associated with the numerical example you wish to reproduce.

---

## Citation

If you use this code or the associated dataset in your work, please cite both the manuscript and the dataset.

### Manuscript

Hijazi, S., Fulgence, N., Burmester, H., Rauter, N., & Gräßle, C.

*Data-driven model order reduction for wave propagation in materials with damage and nonlinearities.*

### Dataset

Hijazi, S., Fulgence, N., Burmester, H., Rauter, N., & Gräßle, C. (2026).

*Data for Data-Driven Model Reduction for Wave Problems* [Data set]. Zenodo.

https://doi.org/10.5281/zenodo.23194087

```bibtex
@dataset{hijazi_wave_problems_2026,
  author       = {Hijazi, S. and Fulgence, N. and Burmester, H. and Rauter, N. and Gr{\"a}{\ss}le, C.},
  title        = {Data for Data-Driven Model Reduction for Wave Problems},
  publisher    = {Zenodo},
  year         = {2026},
  doi          = {10.5281/zenodo.23194087},
  url          = {https://doi.org/10.5281/zenodo.23194087}
}
```

---

## License

This repository is licensed under the **GNU General Public License v3.0 (GPL-3.0)**. See the [`LICENSE`](./LICENSE) file for details.

The associated dataset hosted on Zenodo is made available under the terms of the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**.

---

## Contact

For questions concerning the code or reproduction of the numerical results, please contact:

**Saddam Hijazi**

Institute for Partial Differential Equations  
TU Braunschweig

[saddam.hijazi@tu-braunschweig.de](mailto:saddam.hijazi@tu-braunschweig.de)