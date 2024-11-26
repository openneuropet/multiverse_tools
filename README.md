# Multiverse Tools

Welcome to **Multiverse Tools**, a collection of scripts and resources accompanying the manuscript **"A sensitivity analysis of preprocessing pipelines: Toward a solution for multiverse analyses"**.

This repository provides tools and examples to perform sensitivity analyses on preprocessing pipelines. It includes R code for simulations, a real-world application, and a tutorial on how to carry out a sensitivity analysis.

## Repository Structure

- **`tutorial/`**  
  R code demonstrating how to compute:
  - Global effect estimators
  - Proportion estimator for pipeline sensitivity analysis  
  Refer to **Section 3** of the manuscript for detailed explanations.

- **`sensitivity_analysis-application/`**  
  Code (R and Python) for:
  - Real-world application presented in **Section 4.2**  
  Includes data preprocessing and sensitivity analysis steps.

- **`sensitivity_analysis-simulation/`**  
  R code used for:
  - Simulation works discussed in **Section 4.1**  
  Covers multiple scenarios to test the proposed sensitivity analysis.

## Getting Started

### Installation

Clone this repository:

```bash
git clone https://github.com/openneuropet/multiverse_tools.git
```

Navigate to the desired subdirectory (e.g., `tutorial/`) and follow the instructions in the respective script headers.

### Usage

1. **Run the Tutorial**  
   Explore how to run a sensitivity analysis using simulated data.  
   Navigate to `tutorial/` and run `statistical_sensitivity_analysis_tutorial.Rmd`.

2. **Apply to Real Data**  
   Use the code in `sensitivity_analysis-application/` to replicate the real-world application in **Section 4.2**.

3. **Run Simulations**  
   Generate synthetic data and test the sensitivity analysis using the code in `sensitivity_analysis-simulation/`.

## Contributions

Contributions are welcome! If you identify any issues or have suggestions for improvements, feel free to create a pull request or raise an issue.

## License

This work is licensed under the [MIT License](LICENSE).

## References

- [Ozenne, B; Norgaard, M; Pernet, C; Ganz, M. A sensitivity analysis to quantify the impact of neuroimaging preprocessing strategies on subsequent statistical analyses. 2024.](https://arxiv.org/abs/2404.14882) 

