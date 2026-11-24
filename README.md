# NL_Burger_eq_via_qiskit

One Picard step of the **nonlinear 1D viscous Burgers' equation**, solved with the **HHL quantum linear-systems algorithm** in [Qiskit](https://qiskit.org/).

The nonlinear convection term $u\,u_x$ is handled by **Picard linearization**: the convective velocity is frozen at a previous iterate, which turns one implicit time step into a linear system. That system is non-symmetric, so it is embedded in a **Hermitian dilation** and solved on a 6-qubit HHL circuit simulated with Qiskit Aer. The result is compared against a classical `numpy.linalg.solve` reference.

> **Scope.** This is an educational, small-scale demonstration (4 grid points, 2 clock qubits). The notebook performs a **single** Picard sweep. It shows the linear-solve kernel that HHL would supply inside a full nonlinear (Picard/Newton) loop, but it does not yet iterate to convergence or march in time. It makes no claim of quantum advantage.

A sibling repository with a smaller time step is [Lin_Burger_eq_via_qiskit_HHL](https://github.com/Malik-Saif-Ullah/Lin_Burger_eq_via_qiskit_HHL).

---

## Table of contents

- [Problem statement](#problem-statement)
- [Method](#method)
- [Circuit layout](#circuit-layout)
- [Repository structure](#repository-structure)
- [Getting started](#getting-started)
- [Results](#results)
- [Limitations](#limitations)
- [Possible extensions](#possible-extensions)
- [References](#references)
- [License](#license)

---

## Problem statement

The 1D viscous Burgers' equation on $x \in [0, L]$:

$$
u_t + u\,u_x = \nu\,u_{xx}
$$

| Parameter | Value |
|---|---|
| Domain length $L$ | 1.0 |
| Interior grid points $N$ | 4 ($x = 0.2, 0.4, 0.6, 0.8$), $\Delta x = 0.2$ |
| Viscosity $\nu$ | 0.05 |
| Time step $\Delta t$ | 0.5 (implicit) |
| Boundary conditions | $u(0,t)=1$, $u(L,t)=0$ |
| Previous iterate $u^n$ | linear profile between the BCs plus $0.1\sin(\pi x/L)$ |

The equation is discretized directly (no Cole–Hopf transform).

## Method

### 1. Picard linearization

The convective coefficient is frozen at the previous iterate $u^n$:

$$
u_t + u^n(x)\,u_x = \nu\,u_{xx}
$$

This is linear in the unknown $u$, so one sweep is a single linear solve $A\mathbf{u}^{n+1} = \mathbf{b}$.

### 2. Discretization (backward Euler + central differences)

With the Fourier number $Fo = \nu\Delta t/\Delta x^2$ and the local CFL number $CFL_i = u_i^n \Delta t / (2\Delta x)$:

$$
-(Fo + CFL_i)\,u_{i-1}^{n+1} + (1 + 2Fo)\,u_i^{n+1} - (Fo - CFL_i)\,u_{i+1}^{n+1} = u_i^n
$$

For these parameters $Fo = 0.625$, so the diagonal of $A$ is $2.25$. The first-derivative convection term makes $A$ **non-Hermitian**. Boundary values enter through $\mathbf{b}$.

### 3. Hermitian dilation

HHL needs a Hermitian matrix, so the notebook builds

$$
H = \begin{pmatrix} 0 & A \\ A^\dagger & 0 \end{pmatrix},
\qquad
H \begin{pmatrix} 0 \\ \mathbf{x} \end{pmatrix} = \begin{pmatrix} \mathbf{b} \\ 0 \end{pmatrix}
\;\Longleftrightarrow\; A\mathbf{x} = \mathbf{b}
$$

$H$ is Hermitian for any invertible $A$, and its eigenvalues are $\pm$ the singular values of $A$. The solution sits in the **second half** of the dilated vector, selected by one extra "block" qubit.

### 4. HHL circuit

1. **State preparation:** encode $|\mathbf{b}\rangle$ on the node qubits (block qubit $= |0\rangle$).
2. **Quantum phase estimation (QPE):** controlled $e^{iHt}$ and $e^{2iHt}$ on a 2-qubit clock, then an inverse QFT. The time $t = 0.8\pi/\sigma_{\max}$ keeps the largest $|\lambda|$ below the positive/negative split of the 2-bit signed phase decoding.
3. **Signed eigenvalue inversion:** controlled $R_y$ rotations on the ancilla with angle $\theta = 2\arcsin(C/\lambda_j)$, where clock values are decoded in two's-complement style.
4. **Uncompute:** inverse QPE.
5. **Measure and post-select** on `ancilla = 1` **and** `block = 1`. The node-qubit statistics give $|\mathbf{x}|$.

### 5. Classical comparison

The result is compared with `np.linalg.solve(A, b)`. In the notebook's output this reference column is labelled `u_exact`, but it is the **classical solution of the same linearized system**, not an exact solution of Burgers' equation. HHL returns only a normalized state, and amplitudes are recovered as $\sqrt{p}$ from measured probabilities, so the notebook rescales them by the norm of the classical solution.

## Circuit layout

| Register | Qubits | Role |
|---|---|---|
| `sys[0]`, `sys[1]` | 2 | Node index (4 grid points) |
| `sys[2]` | 1 | Dilation block qubit |
| `clock` | 2 | Phase-estimation register |
| `ancilla` | 1 | Eigenvalue-inversion flag |
| **Total** | **6** | |

Classical register (4 bits): `ancilla`, `block`, `node0`, `node1`.

## Repository structure

```
NL_Burger_eq_via_qiskit/
├── NL_Burger_eq_via_qiskit.ipynb   # Full pipeline, with executed outputs
├── requirements.txt                # Python dependencies
├── LICENSE                         # MIT
├── .gitignore
└── README.md
```

## Getting started

**Requirements:** Python 3.10+ (developed on 3.11) and Qiskit ≥ 2.1.

```bash
git clone https://github.com/Malik-Saif-Ullah/NL_Burger_eq_via_qiskit.git
cd NL_Burger_eq_via_qiskit

python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt

jupyter notebook NL_Burger_eq_via_qiskit.ipynb
```

Run the single code cell. It prints each step (setup, matrix, dilation, eigenvalues, decoding, counts, comparison).

## Results

Example run (Qiskit 2.6.2, Aer 0.18.4, 200,000 shots). **Post-selection success probability: ≈ 61.9%.**

Spectrum of $H$: $\pm\{1.560,\ 2.203,\ 2.802,\ 3.170\}$.

| Node | x | u (HHL) | u (classical) | Abs. error | Rel. error |
|---|---|---|---|---|---|
| 1 | 0.2 | 1.0254 | 0.9628 | 0.0626 | 6.5% |
| 2 | 0.4 | 0.9350 | 0.8719 | 0.0631 | 7.2% |
| 3 | 0.6 | 0.5920 | 0.7032 | 0.1112 | 15.8% |
| 4 | 0.8 | 0.2735 | 0.4114 | 0.1379 | 33.5% |

HHL captures the overall decaying profile, but the error grows toward the outflow boundary and reaches about 34% at the last node. The measurement is stochastic and no random seed is set, so values change slightly between runs.

Compared with the sibling repository ($\Delta t = 0.1$, errors of roughly 4–7%), the larger time step here gives a stronger-coupled, wider-spectrum system and noticeably larger errors. The likely reason is the coarse phase estimation described below, though this has not been isolated experimentally.

## Limitations

- **One Picard sweep.** The nonlinearity is only linearized once about a fixed guess. The notebook does not iterate, check convergence, or advance in time.
- **Coarse phase estimation.** The dilated matrix has 8 distinct eigenvalues (magnitudes 1.56 to 3.17), but a 2-qubit clock has only 4 phase bins, of which 3 are used for inversion. Eigenvalues are therefore not resolved individually.
- **Hamiltonian simulation by matrix exponentiation.** $e^{iHt}$ is built classically with `scipy.linalg.expm` and wrapped as a unitary. This is fine for a demo but does not scale and does not reflect real HHL gate costs.
- **Magnitudes only.** Amplitudes are $\sqrt{p}$ from measurement statistics, so signs are lost. The solution here is positive.
- **Classical rescaling.** The overall scale of the HHL output uses the classical solution's norm. Reading out the full vector also removes the usual exponential-speedup argument for HHL.
- **Simulator only.** No noise model or real hardware.

## Possible extensions

- Wrap the solver in a full Picard loop (update $u^n \leftarrow u^{n+1}$ until convergence) and time-march toward a shock
- Increase $N$ and the number of clock qubits and study error versus resource scaling
- Replace the exact unitary with Trotterized or block-encoded Hamiltonian simulation
- Add a noise model or run on real hardware
- Compare against variational quantum linear solvers

## References

1. A. W. Harrow, A. Hassidim, S. Lloyd, "Quantum algorithm for linear systems of equations," *Physical Review Letters* **103**, 150502 (2009).
2. J. M. Burgers, "A mathematical model illustrating the theory of turbulence," *Advances in Applied Mechanics* **1**, 171–199 (1948).
3. [Qiskit documentation](https://docs.quantum.ibm.com/) and [Qiskit Aer](https://qiskit.github.io/qiskit-aer/).

## License

Released under the [MIT License](LICENSE).
