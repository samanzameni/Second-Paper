# Response-Level Differential Privacy for Geo-Location Databases: Power Cost, Coalition Limits and Comparison with Input Perturbation

**Saman Zameni, MohammadReza Taheri, Hadi Pouladi**
**Status:** Submitted to *Ad Hoc Networks*

> [!IMPORTANT]
> **This is not the final code for this paper.** The repository contains a preliminary version of the simulation code. The final code used to produce the results reported in the paper will be released here once the paper is published.

## Abstract

Maximum transmission power (MTP) limits returned by geo-location databases (GDBs) depend on, and can reveal, the positions of protected primary users (PUs). We study DP-MTP, a response-level mechanism that adds Laplace noise calibrated to a Lipschitz bound of the MTP ceiling, subtracts a public backoff and clips. It is *ε*-differentially private for PU displacements of up to *R* metres and bounds the violation probability of every fresh release; before clipping, its power loss is at most 8.5% above a lower bound for additive, translation-invariant metric-private releases at a 1% violation target. Van Trees bounds show coalition accuracy grows with the square root of the number of answers; a public per-cell ledger with cached releases bounds each PU’s loss per epoch for any number of colluding devices. In synthetic experiments (*R* = 1 m, *ε* = 1), one response raises the distance error of informative queries from zero to 13.0 m at 91.5% retained power, yet adaptive coalitions localize a PU to 0.28 m after 200 queries and defeat cloaking, sybil-assisted *k*-anonymity and coarsening. With the ledger, tested attacks stay near 2 m within one epoch but reach 0.43 m after eight. Perturbing each registered PU position once per epoch instead certifies an error of at least 2.0 m for every coalition (2.4 m for the Bayes-optimal one) with a per-epoch budget of 1 instead of 45, retaining 95% of the power. For static incumbents, registration-level perturbation is preferable; response-level noise may pay off when few responses depend on a PU between registration updates, a hypothesis left for testing.

**Keywords:** Differential privacy, Geo-location database, Cognitive radio, Spectrum sharing, Location privacy, Metric privacy

## Repository contents

| File | Purpose |
|---|---|
| `DP_MTP_AllFigures.m` | Main script: sets the parameters, runs all experiments and saves the figures |
| `run_mc.m` | Monte Carlo simulation of the DP-MTP mechanism (noise error, retained power, attack RMSE) |
| `run_baselines.m` | Baseline privacy mechanisms simulated for comparison with DP-MTP |
| `particle_attack.m` | Particle/importance-sampling attack used to estimate PU distance from noisy responses |
| `Fig*.png` | Figures produced by the main script |

## How to run

Requires MATLAB (R2021a or later) with the Statistics and Machine Learning Toolbox.

```matlab
DP_MTP_AllFigures   % runs run_mc, run_baselines and particle_attack, then saves Fig1–Fig7 (rng seed 42)
```

## Example outputs (preliminary code)

| Attack RMSE vs. ε | Comparison with baselines |
|---|---|
| ![RMSE vs epsilon](Fig6_RMSE_vs_epsilon.png) | ![Baseline comparison](Fig7_Baseline_Comparison_Bar.png) |

These figures come from this preliminary version and may differ from the figures in the paper.

## Contact

Saman Zameni · [Website](https://samanzameni.netlify.app) · [ORCID](https://orcid.org/0009-0006-8774-6391) · [saman.zameni@gmail.com](mailto:saman.zameni@gmail.com)
