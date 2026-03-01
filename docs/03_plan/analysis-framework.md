# Analysis Framework — Phase 3

## Research Question (from Phase 2)

In the exactly solvable 1D TFIM (and anisotropic XY chain), does a quantum quench from an initial
product state |ψ_θ⟩ (that breaks Z2 symmetry at angle θ) to H(g_f) exhibit the Z2 Quantum Mpemba
Effect (EA relaxation curve of more-asymmetric state crossing below less-asymmetric state at finite
time t_M), and does this occurrence **coincide with the quench crossing the DQPT boundary**?

Three sub-questions locked in Phase 2:
1. (a) QME occurs ↔ quench crosses DQPT boundary (g_i vs. g_f straddle g_c)
2. (b) Mpemba crossing time t_M diverges critically as g_f → g_c
3. (c) QME ↔ DQPT correspondence is separable from QPT using the anisotropic XY chain

---

## Models

### Model 1: Transverse-Field Ising Model (TFIM)
H = -J Σ_i σ^z_i σ^z_{i+1} - g Σ_i σ^x_i   (J=1)

Via Jordan-Wigner, maps to BdG free-fermion Hamiltonian.
- QPT boundary: g_c = 1
- In TFIM, QPT and DQPT coincide (g_i < 1, g_f > 1 ↔ DQPT occurs)
- Problem: QPT and DQPT boundaries are not independently tuneable in TFIM

### Model 2: Anisotropic XY Chain (DECOUPLING model)
H = -Σ_i [(1+γ)/4 σ^x_i σ^x_{i+1} + (1-γ)/4 σ^y_i σ^y_{i+1}] - h/2 Σ_i σ^z_i

Also maps to BdG free-fermion via Jordan-Wigner.
- QPT boundary: h_c = 1 (at fixed anisotropy γ)
- DQPT structure: Depends on (h_i, h_f, γ) — CAN differ from QPT boundary
- Key: For certain (γ, h_i, h_f), DQPT occurs even for quenches within the same QPT phase
  → Allows DIRECT TEST of QME ↔ DQPT hypothesis independent of QPT

---

## Initial States

Family of Z2-symmetry-breaking product states parameterized by angle θ:
|ψ_θ⟩ = Π_i [cos(θ/2)|↑⟩_i + sin(θ/2)|↓⟩_i]

- θ = 0: fully polarized (maximal Z2 breaking)
- θ = π/2: equal superposition
- Different θ → different initial EA_0(θ)
- Two states θ₁ < θ₂ with EA_0(θ₁) > EA_0(θ₂) (or vice versa) are chosen for QME test

Under Jordan-Wigner, |ψ_θ⟩ maps to a specific Gaussian fermionic state with
known covariance matrix Γ(0).

---

## Observables

### 1. Entanglement Asymmetry (EA) — Primary Observable
ΔS_A(t) = S(ρ̃_A(t)) - S(ρ_A(t))

where ρ̃_A = (ρ_A + Z_A ρ_A Z_A†) / 2  (Z2 twirl, Z_A = Π_{j∈A} σ^z_j)

**Computation via Majorana covariance matrix:**
- Full covariance matrix Γ(0) for initial state (analytically derived)
- Time evolution: Γ(t) = e^{K t} Γ(0) e^{-K t}  (K = BdG generator, skew-symmetric)
- Restrict to subsystem A: Γ_A(t) = Γ(t)|_{A×A}
- S(ρ_A): from eigenvalues ν_k of iΓ_A → S = Σ_k h((1+ν_k)/2)
- S(ρ̃_A): use parity-projected Gaussian state formula
  - ⟨Z_A⟩ = Pf(Γ_A) or via det(1-2C_A) where C_A is complex correlator submatrix
  - p_± = (1 ± ⟨Z_A⟩) / 2
  - S(ρ̃_A) = H(p_+, p_-) + p_+ S_+  + p_- S_-  (symmetry-resolved entropies of ±sectors)
  - S_± from modified boundary-condition correlators (standard BdG technique)

### 2. Loschmidt Echo — DQPT Diagnostic
L(t) = |⟨ψ_0|e^{-iH_f t}|ψ_0⟩|²
g(t) = (1/N) log L(t)  (return rate function)

DQPTs identified as non-analytic cusps in g(t) at times t*.

**Computation:** L(t) = |det(cos(ε_k t) + i sin(ε_k t) cos(Δ_k))|² (product over k)
where ε_k = single-particle energies, Δ_k = Bogoliubov angle mismatch.

### 3. Mpemba Crossing Time t_M
t_M(θ₁, θ₂; g_f) = min{t > 0 : ΔS_A(θ₁, t) = ΔS_A(θ₂, t)}
(with ΔS_A(θ₁, 0) > ΔS_A(θ₂, 0))

t_M → ∞ means no QME for this parameter combination.

---

## Analysis Dimensions and Story Lines

### Main Story Line (Paper Argument)
> DQPTs act as a control switch for the quantum Mpemba effect: QME occurs *if and only if*
> the quench dynamics generates DQPTs. This is demonstrated exactly in free-fermion models,
> revealed via the Majorana covariance matrix calculation of Z2 entanglement asymmetry.

### Supporting Lines

**SL1 — QME Phase Diagram in TFIM**
- Scan (g_i, g_f) ∈ [0.2, 3]² on a 10×10 grid; for each pair, run two θ values
- Determine: QME present or absent for each (g_i, g_f)
- Overlay DQPT boundary onto QME phase diagram
- Expected: QME region = DQPT region

**SL2 — EA Dynamics and Crossing**
- Show representative EA(t) curves for (g_i, g_f) in:
  (a) same phase, no DQPT → no crossing (no QME)
  (b) different phases, DQPT → crossing (QME)
- Show Loschmidt echo g(t) alongside to demonstrate DQPT timing
- Expected: first EA crossing correlates with first DQPT time t*

**SL3 — Critical Behavior of t_M**
- Fix g_i < 1, scan g_f from 1.0 → 3.0
- Plot t_M(g_f): expected divergence as g_f → 1+ (QC boundary)
- Fit: t_M ~ (g_f - 1)^{-α} → extract critical exponent α
- Also scan from below (g_i > 1 → g_f < 1) for symmetry
- Expected: t_M diverges at DQPT boundary, no QME at DQPT-free quenches

**SL4 — DQPT/QPT Decoupling via XY Chain (Key Test)**
- Use XY chain with γ ∈ {0.3, 0.5, 0.7}
- Find parameter region where DQPT occurs within SAME QPT phase (quench within ordered phase, but crosses DQPT line)
- Run QME test: does QME occur for within-phase quench with DQPT?
- Expected: YES — QME follows DQPT, NOT QPT → proves the DQPT mechanism

**SL5 — Finite-Size Scaling (Robustness)**
- Run main results at N ∈ {100, 500, 1000, 2000, 5000}
- Show QME signal persists / strengthens as N → ∞
- Finite-size scaling of t_M near critical boundary
- Expected: QME is a genuine thermodynamic phenomenon

**SL6 — Analytical Physical Argument (quasiparticle picture; framed as conjecture + evidence)**
- Physical mechanism: DQPT-crossing quenches cause quasiparticle mode occupations n_k(t)
  to cross 1/2 at DQPT critical times t*. This causes parity ⟨Q_A(t)⟩ = ∏_k(1-2n_k) to
  oscillate in sign. For non-DQPT quenches, all n_k stay bounded away from 1/2 → no sign flip.
- Conjectured condition (to verify): QME occurs ↔ ∃k ∈ BZ: Δ_k > π/2 (DQPT topological condition)
- Frame: "We derive an analytical physical argument based on quasiparticle picture and verify numerically"
  NOT: "We prove an equivalence theorem" (would require more rigorous proof)
- Literature context:
  * Liu et al. (PRL 2024, arXiv:2403.08459): Z2 QME absent in random (chaotic) circuits
  * This paper: Z2 QME occurs in integrable systems specifically when DQPT is present
  * DQPT supplies coherent quasiparticle dynamics absent in chaotic systems → enables Z2 QME
  * This contrast between integrable and non-integrable gives physical interpretation

**SL7 — Control: Ordinary Entanglement Entropy (new, per Domain Scientist)**
- Compute S(ρ_A)(t) for same quenches as SL2 — confirm NO crossing (no generic Mpemba)
- Shows the QME is SPECIFIC to Z2 symmetry sector, not a generic feature of entanglement
- One figure panel; simple to compute

---

## Computation Plan (revised after Round-1 agent review)

### Software Stack
- Python 3.12, NumPy, SciPy (exact free-fermion computation)
- joblib.Parallel for 128-CPU parallelism (no GPU needed — 2-5 min total runtime)
- Matplotlib for figures

### Time Evolution — Implementation Rules (after Methodology Consultant review)
- BdG diagonalization: O(N³) once per Hamiltonian. Store eigenvectors V and eigenvalues ε_k.
- Single-step time evolution: C(t) = V · diag(e^{-iε_k t}) · (V† C(0) V) · diag(e^{iε_k t}) · V†
  NO step-by-step matrix multiplication (accumulates orthogonality drift).
- Subsystem C_A(t): extract N_A × N_A block from C(t) after full evolution.
- det(1-2C_A): use eigenvalues of C_A, compute as product ∏_k(1-2λ_k), never np.linalg.det.

### S± Computation (corrected after Round-2 Methodology Consultant review)
CORRECTED FORMULA (two bugs fixed from Round-1 proposal):

Bug 1 fixed: r(1/2) = 1 - ln2 ≠ 0 (l'Hôpital applied correctly below)
Bug 2 fixed: S_± formula includes + log(p_±) term

Define per-mode functions:
  a_k = λ_k log λ_k + (1-λ_k) log(1-λ_k)                 [entropy weight; ≤ 0]
  r_k = [(1-λ_k) log(1-λ_k) - λ_k log λ_k] / (1-2λ_k)    [parity-weighted]
       = (2ln2 - 2)/(-2) = 1-ln2  when λ_k = 1/2  (l'Hôpital: f'(λ)/g'(λ) at λ=1/2)
       where f'(λ) = -ln(1-λ) - ln(λ) - 2,  g'(λ) = -2
       → r(1/2) = (2ln2 - 2)/(-2) = 1 - ln2 ≈ 0.307

Partial derivatives at n=1:
  ∂_n P₀|_{n=1} = Σ_j a_j  (= -S(ρ_A))
  ∂_n P_π|_{n=1} = P_π · Σ_j r_j  where P_π = ⟨Q_A⟩ = ∏_k(1-2λ_k)

CORRECT S± formula (includes + log p_± for normalization):
  S_+ = -[∂_n P₀ + ∂_n P_π] / (2 p_+)  +  log p_+
  S_- = -[∂_n P₀ - ∂_n P_π] / (2 p_-)  +  log p_-

Guard at p_- → 0 (near DQPT critical time t*): set S_- = 0 explicitly.

Sanity checks:
  At θ=π/2 (λ_k = 1/2 for all k): P_π=0, p_±=1/2, S_± = S(ρ_A) - ln2  ✓ (equipartition)
  At θ=0 (λ_k=0 for all k): P_π=1, p_+=1, S_+=S(ρ_A), ΔS_A=0  ✓ (Z2-symmetric state)

Validated against small-N_A exact diagonalization (N_A = 4, 6, 8).

### Revised Time Window (corrected after Round-2 Statistical Rigor review)
N-dependent t_max (required for finite-size scaling runs with N=5000):
  t_max(N) = max(60/J,  3*N_A / (4*J))   [Lieb-Robinson saturation criterion]
  N=100:  t_max = max(60, 7.5) = 60/J
  N=500:  t_max = max(60, 37.5) = 60/J
  N=1000: t_max = max(60, 75) = 75/J
  N=2000: t_max = max(60, 150) = 150/J
  N=5000: t_max = max(60, 375) = 375/J

Time resolution:
  dt = 0.02 for t ∈ [0, 20/J]
  dt = 0.05 for t ∈ [20/J, t_max]

QME classification threshold: RELATIVE 1% rule (replacing absolute 10^{-4}):
  QME PRESENT if: min_t [ΔS_A(θ₁,t) - ΔS_A(θ₂,t)] / |ΔS_A(θ₁,0) - ΔS_A(θ₂,0)| < -0.01
  Borderline if reversal < 1% of initial gap → flagged for extended time window check
  Pairs with |ΔS_A(θ₁,0) - ΔS_A(θ₂,0)| < 0.05 excluded from universality count

Verify EA is within 1% of saturation value at t_max for all (N, g_i, g_f) combinations.

### Full Parameter Scan (revised)
- TFIM QME phase diagram: 10×10 coarse grid + 20×20 fine grid near g=1 (g∈[0.5,1.5])
  × 4 θ pairs = ~1200 runs (~5 min with 128-CPU parallelism)
- θ-pair independence: 5×5 grid in (θ₁,θ₂) ∈ [π/8, 7π/8]² at 4 representative (g_i,g_f) points = 100 runs
- Critical exponent scan: 1 g_i × 100 g_f (dense near g=1) × 3 θ = 300 runs
- XY chain: 3 γ × 8×8 grid × 4 θ pairs = 768 runs
- Finite-size scaling: N ∈ {100, 200, 500, 1000, 2000, 5000} × 6 representative points = 36 runs
- N_A/N sensitivity: N_A/N ∈ {0.05, 0.1, 0.15, 0.2} at 4 parameter points = 16 runs
Total estimate: ~15 min on 128 CPUs

### Validation Checks (expanded after Round-1 review)
1. ΔS_A(t=0, θ, N_A): closed-form formula for product state (all λ_k = sin²(θ/2))
   - ⟨Z_A⟩ = cos^{N_A}(θ);  p± = (1 ± cos^{N_A}θ)/2;  S(ρ_A) = N_A h(sin²(θ/2))
2. ΔS_A = 0 for θ = 0 (Z2-symmetric initial state — exact zero)
3. S(ρ_A) at g=1 (CFT): S = (1/6) log N_A + const  (c=1/2)
4. Loschmidt echo: compare t* to analytical TFIM formula (Fisher zeros at k*)
5. Loschmidt phase winding number: track arg G(t) for topological DQPT signature
6. Rényi EA ΔS_A^{(2)} (n=2): compare to charged-moment formula — cross-check S±
7. Small N_A exact diagonalization: N_A=4,6,8 — compare every formula to explicit matrix calc
8. Thermodynamic limit: Toeplitz matrix + Szegő theorem for N→∞ check

### Additional Observables (added after Round-1 review)
- Ordinary entanglement entropy S(ρ_A) alongside ΔS_A — confirm NO Mpemba crossing for S_A
  (isolates Z2-specific nature of QME)
- Loschmidt echo phase arg G(t) and winding number (topological DQPT indicator)
- Entanglement spectrum {λ_k(t)} vs. time — dynamical signatures
- Site-resolved magnetization ⟨σ^z_i(t)⟩ for physical interpretation

### QME Classification Rule (explicit, per Statistical Rigor Advisor)
QME is declared PRESENT for parameter pair (g_i, g_f, θ₁, θ₂) if:
  min_{t ∈ (0, t_max]} [ΔS_A(θ₁, t) - ΔS_A(θ₂, t)] < 0
  where ΔS_A(θ₁, 0) > ΔS_A(θ₂, 0)  (θ₁ is more asymmetric initially)
QME is ABSENT if: ΔS_A(θ₁, t) > ΔS_A(θ₂, t) for all t ∈ (0, t_max]
Borderline cases (relative reversal < 1% of initial EA gap): flagged for extended time window.
Pairs with initial EA gap < 0.05 excluded from phase diagram classification.

Pre-registered θ-grid outcome criteria (from Statistical Rigor Advisor):
  STRONG CONFIRMATION: QME present for ≥90% of θ pairs in DQPT-crossing quenches;
                        ABSENT for ≥90% in non-DQPT-crossing quenches → "iff" claim holds
  WEAK CONFIRMATION: QME present for 60-89% → qualify claim as "generically correlated"
  FALSIFICATION: QME present for any θ pair in non-DQPT quench → "iff" must be dropped

---

## Sufficiency Criteria (D-Step 2)

- [ ] QME phase diagram in TFIM (10×10 grid) computed
- [ ] Both DQPT-crossing and DQPT-absent quenches demonstrated
- [ ] t_M divergence extracted with critical exponent
- [ ] XY chain test showing QME ↔ DQPT (not QPT)
- [ ] Finite-size scaling from N=100 to N=5000
- [ ] At least one analytical result connecting QME condition to DQPT topology
- [ ] 4+ main figures, 1+ table, 20+ references
- [ ] All results validated against known benchmarks

---

## Alternative Explanations (D-Step 3)

| Expected Finding | Alternative Explanation | How to Rule Out |
|---|---|---|
| QME occurs ↔ DQPT | QME depends on initial EA magnitude only (larger θ always wins) | Show: same θ₁,θ₂ pair can have QME or no-QME depending only on g_f crossing the DQPT boundary |
| QME ↔ DQPT, not QPT | QME tracks QPT (equilibrium phase transition), not the dynamical boundary | XY chain SL4: find DQPT within same QPT phase → show QME still occurs |
| t_M diverges at DQPT boundary | t_M divergence is finite-size artifact | Finite-size scaling SL5: show divergence sharpens with N |
| t_M ~ (g_f - g_c)^{-α} power law | t_M has logarithmic divergence (not power law) | Fit both power law and log law; compare BIC/AIC; verify with analytical prediction |
| QME is a universal phenomenon tied to DQPT universality class | QME is model-specific (TFIM only) | Verify in XY chain with different universality class (Ising vs. XY) |

---

## Story Line Narrative

**Title (candidate):** "Quantum Mpemba Effect Controlled by Dynamical Quantum Phase Transitions"

**Abstract arc:**
1. Motivation: The quantum Mpemba effect — counterintuitive accelerated symmetry restoration — lacks
   a universal control mechanism. When does it occur?
2. Question: We ask whether dynamical quantum phase transitions (DQPTs) in 1D integrable systems
   act as the control mechanism.
3. Method: We use the Majorana covariance matrix to compute exactly the Z2 entanglement asymmetry
   in the TFIM and anisotropic XY chain after quantum quenches.
4. Result: QME occurs if and only if the quench generates DQPTs. The Mpemba crossing time t_M
   diverges critically at the DQPT boundary with a universal exponent. This QME↔DQPT link is
   distinct from the QPT boundary.
5. Significance: Establishes DQPTs as a dynamical control parameter for quantum information
   scrambling/restoration, with implications for quantum thermalization.
bling/restoration, with implications for quantum thermalization.
