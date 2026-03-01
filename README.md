# Z2 Quantum Mpemba Effect in Integrable Free-Fermion Chains

> **Decoupling from Dynamical Quantum Phase Transitions**

---

## ⚠️ Important Notice

> This research project was completed **entirely** by [**Amplify**](https://evoclaw.github.io/amplify/) —
> an agentic research automation framework — with **no human editing** of scientific content,
> code, or manuscript.
>
> The results are **reproducible** from the provided code and data.
> Nevertheless, as an AI-generated research output, **please treat all findings with appropriate caution**
> and verify independently before building upon them.

---

## 📄 Paper

### **[→ Download / View PDF](paper/main.pdf)**

> *Z2 Quantum Mpemba Effect in Integrable Free-Fermion Chains: Decoupling from Dynamical Quantum Phase Transitions*
>
> **Author:** [Amplify](https://evoclaw.github.io/amplify/) (agentic research framework)

---

## Overview

This project investigates the **Z2 quantum Mpemba effect (QME)** in one-dimensional integrable
free-fermion models — the transverse-field Ising model (TFIM) and the anisotropic XY chain.

The quantum Mpemba effect is the counter-intuitive phenomenon where a *more* strongly
symmetry-broken quantum state restores its symmetry *faster* than a less broken companion.
We study the discrete Z2 case using the **entanglement asymmetry** as a probe.

### Key Findings

| Finding | Result |
|---|---|
| Z2 QME in integrable systems | **Universal** — present at all 23 tested post-quench fields, including the DQPT-free FM phase |
| Relation to DQPTs | **Decoupled** — DQPTs modulate but do not govern crossing times |
| Predictor of crossing time $t_M$ | Initial EA imbalance $\|\Delta\Delta\mathcal{S}\|$ (Spearman $\rho = +0.90$) |
| Thermodynamic limit | $t_M$ converges to finite values (FSS over $N = 8$–$18$) |
| Model robustness | Present across full XY chain family ($\gamma = 0$–$1$) |

### Exact Formula for Initial Entanglement Asymmetry

$$\Delta S_A(0;\theta) = h_{\rm bin}\!\left(\frac{1 + \cos^{N_A}\theta}{2}\right)$$

where $h_{\rm bin}(p) = -p\ln p - (1-p)\ln(1-p)$.  
Verified numerically to relative error $< 10^{-14}$.

---

## Repository Structure

```
z2-quantum-mpemba/
├── paper/
│   ├── main.pdf              ← compiled paper (start here)
│   ├── main.tex              ← master LaTeX file
│   ├── preamble.tex
│   ├── references.bib        ← 18 verified references
│   ├── sections/             ← abstract, intro, methods, results, discussion, conclusion
│   └── supplementary/
├── src/
│   ├── exp_common.py         ← core physics library (TFIM, EA, DQPT detection)
│   ├── exp_a_v2.py           ← Exp A: dynamics in FM and PM phases
│   ├── exp_b_v2.py           ← Exp B: 671-pair statistical analysis
│   ├── exp_c_v2.py           ← Exp C: finite-size scaling
│   ├── supp_r1_larger_N.py   ← R1: FSS up to N=18
│   ├── supp_r2_gf_scan.py    ← R2: 23-value g_f scan
│   ├── supp_r4_xy_chain.py   ← R4: anisotropic XY chain
│   ├── supp_s2_statistics.py ← S2: KS test and binomial proximity test
│   ├── plot_all_v2.py        ← generate all main figures
│   └── plot_supplements.py   ← generate supplementary figures
├── results/
│   ├── figs/                 ← all publication figures (PDF + PNG)
│   ├── exp_b_v2.json         ← main dataset: 671 Mpemba pairs
│   ├── supp_r1_fss.json      ← FSS results (N=8–18)
│   ├── supp_r2_gf_scan.json  ← g_f scan (23 values)
│   ├── supp_r4_xy.json       ← XY chain results
│   └── supp_s2_stats.json    ← statistical test results
└── docs/
    └── 06_integration/
        └── argument-blueprint.md  ← paper argument structure
```

---

## Reproducing the Results

### Requirements

```bash
pip install numpy scipy matplotlib
```

Python 3.9+ recommended. No GPU required; all computations run on CPU.

### Run experiments (in order)

```bash
# Core experiments
python src/exp_a_v2.py          # FM vs PM dynamics  (~2 min, N=12)
python src/exp_b_v2.py          # 671-pair statistics (~5 min, N=12)
python src/exp_c_v2.py          # finite-size scaling (~3 min, N=8–14)

# Supplementary
python src/supp_r1_larger_N.py  # FSS up to N=18    (~10 min)
python src/supp_r2_gf_scan.py   # 23-field scan      (~5 min, N=14)
python src/supp_r4_xy_chain.py  # XY chain           (~5 min)
python src/supp_s2_statistics.py # statistical tests  (<1 min)

# Generate all figures
python src/plot_all_v2.py
python src/plot_supplements.py
```

Results are saved to `results/` as JSON; figures to `results/figs/`.

### Compile the paper

```bash
cd paper
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

---

## Figures

| Figure | File | Description |
|---|---|---|
| Fig. 1 | `fig0_typical_DS` | Representative $\Delta S_A(t)$ with QME crossing |
| Fig. 2 | `fig3_gf_compare` | $t_M$ vs post-quench field $g_f$ (FM and PM phases) |
| Fig. 3 | `fig2_tM_vs_dDS` | $t_M$ vs initial EA imbalance (671 pairs) |
| Fig. 4 | `fig3_fractional_test` | Fractional crossing position $\phi$ distribution |
| Fig. 5 | `fig4_finite_size` | Finite-size scaling of $t_M$ |

---

## About Amplify

This project was produced by **[Amplify](https://evoclaw.github.io/amplify/)**, a structured
agentic research automation framework that guides an AI agent through the full scientific
workflow: literature review → problem validation → method design → experiment execution →
results integration → paper writing.

Amplify enforces scientific rigour through phase-gated progression, multi-agent deliberation,
metric locking, anti-cherry-picking, and reproducibility-driven experiment management.

> This repository serves as a demonstration of Amplify's capabilities.
> All code was written, all experiments were run, and all manuscript text was drafted
> autonomously by the AI agent under Amplify's workflow.
> **No human scientific editing was performed.**

---

## Citation

If you use this work, please cite:

```bibtex
@misc{z2qme2025,
  title  = {Z2 Quantum {Mpemba} Effect in Integrable Free-Fermion Chains:
            Decoupling from Dynamical Quantum Phase Transitions},
  author = {Amplify},
  year   = {2025},
  note   = {AI-generated research via Amplify agentic framework.
            \url{https://evoclaw.github.io/amplify/}},
  url    = {https://github.com/EvoClaw/z2-quantum-mpemba}
}
```

---

## License

Code: [MIT License](LICENSE)

*Paper and figures: generated by Amplify for demonstration purposes.
Treat with caution — independent verification is recommended before scientific use.*
