# Mathematical and Biomechanical Foundations of Local Dynamic Stability (LDS) and Largest Lyapunov Exponent (LLE) Estimation

## 1. Introduction and Dynamical Systems Formulation

Human locomotion is governed by high-dimensional, nonlinear, dissipative neuromusculoskeletal dynamics. While traditionally analyzed via discrete spatiotemporal metrics (e.g. stride length variability), the underlying stability is best characterized through **dynamical state-space reconstruction** and **finite-time Lyapunov exponents**.

A continuous dynamical system on a manifold $\mathcal{M}$ is defined by:
$$\dot{\mathbf{x}}(t) = \mathbf{F}(\mathbf{x}(t)), \quad \mathbf{x} \in \mathbb{R}^D$$

For an autonomous attractor $\mathcal{A}$, the Largest Lyapunov Exponent $\lambda_1$ measures the mean exponential rate of divergence of infinitesimally close trajectories:
$$d(t) \approx d(0) e^{\lambda_1 t} \implies \lambda_1 = \lim_{t \to \infty} \lim_{d(0) \to 0} \frac{1}{t} \ln \frac{d(t)}{d(0)}$$

In experimental biomechanics, the true state vector $\mathbf{x}(t)$ is unobservable; only a scalar observation $s(t) = h(\mathbf{x}(t))$ is recorded (e.g. joint angle $\theta(t)$ or marker velocity $v(t)$).

---

## 2. Takens' Delay Embedding Theorem (1981)

Takens' theorem guarantees that for a smooth generic dynamical system on a compact manifold of dimension $d_A$, the delay coordinate map:
$$\Phi_{(\tau, d_E)}(s(t)) = \mathbf{S}(t) = \left[ s(t), s(t+\tau), s(t+2\tau), \dots, s(t+(d_E-1)\tau) \right]^T \in \mathbb{R}^{d_E}$$
is an **embedding** (a smooth diffeomorphism onto its image) provided that:
$$d_E > 2 d_A$$

### 2.1. Selection of Delay ($\tau$)
The delay $\tau$ must balance two opposing criteria:
- If $\tau$ is too small: $s(t) \approx s(t+\tau)$ (redundancy, trajectories collapse to the main diagonal).
- If $\tau$ is too large: nonlinear stretching decorrelates components completely (irrelevance).

Following Fraser & Swinney (1986), $\tau$ is chosen as the **first local minimum of the Average Mutual Information (AMI)** function:
$$I(\tau) = \sum_{i,j} P(s_i, s_{i+\tau}) \log_2 \frac{P(s_i, s_{i+\tau})}{P(s_i) P(s_{i+\tau})}$$
In human running normalized to 100 samples/stride, $\tau \approx 10$ samples.

### 2.2. Selection of Embedding Dimension ($d_E$)
Following Kennel, Brown, & Abarbanel (1992), $d_E$ is determined using the **Global False Nearest Neighbors (GFNN)** metric, identifying when projected overlapping orbits uncross in higher dimensions:
$$\frac{|s(t+d_E\tau) - s(t_j+d_E\tau)|}{\|\mathbf{S}_{d_E}(t) - \mathbf{S}_{d_E}(t_j)\|} > R_{tol}$$
In human running kinematics, $d_E \approx 5$ is canonically observed.

---

## 3. Rosenstein's Algorithm and Multi-Neighbor Revision

### 3.1. Classic Formulation (Rosenstein et al., 1993)
For each reconstructed state vector $\mathbf{S}_i \in \mathbb{R}^{d_E}$:
1. Search for its single nearest neighbor $\mathbf{S}_{\hat{j}}$:
   $$\hat{j} = \arg\min_j \|\mathbf{S}_i - \mathbf{S}_j\|, \quad \text{subject to } |i - j| > W$$
   where $W$ is the **Theiler temporal exclusion window**.
2. Track Euclidean divergence forward in discrete steps $k = 1, 2, \dots, K$:
   $$d_i(k) = \|\mathbf{S}_{i+k} - \mathbf{S}_{\hat{j}+k}\|$$
3. Average the logarithm of divergence across all reference points $N$:
   $$y(k) = \langle \ln d_i(k) \rangle = \frac{1}{N} \sum_{i=1}^N \ln d_i(k)$$
4. The slope of $y(k)$ over the linear scaling interval defines $\lambda_1$:
   $$y(k) \approx \lambda_1 \cdot k \cdot \Delta t + C$$

### 3.2. Multi-Neighbor Extension (Mehdizadeh, 2019)
Tracking only $k=1$ nearest neighbor makes the estimate exquisitely fragile to high-frequency measurement noise. Mehdizadeh (2019) extended the formulation to average across the $k_{nn}$ nearest neighbors:
$$y(k) = \frac{1}{N \cdot k_{nn}} \sum_{i=1}^N \sum_{m=1}^{k_{nn}} \ln \|\mathbf{S}_{i+k} - \mathbf{S}_{j(i, m)+k}\|$$
For $k_{nn} \ge 15$, the estimate reaches an asymptote/plateau, neutralizing additive Gaussian and measurement noise.

---

## 4. The Theiler Window vs. Embedding Delay Fallacy

A prevalent error in historical software (e.g. `LyE_R.m` in the NONAN toolbox) is defining:
$$W = 0.8 \tau \quad (\text{INCORRECT})$$

### Why this is mathematically invalid:
Theiler (1986) demonstrated that if points are separated by less than the autocorrelation time $t_{corr}$, they belong to the **same continuous trajectory segment**. Their short-term separation reflects orbital velocity along the flow $\mathbf{F}(\mathbf{x})$, not Lyapunov exponent divergence:
$$\mathbf{S}_{i+\delta} - \mathbf{S}_i \approx \delta \cdot \dot{\mathbf{S}}_i$$
In running kinematics sampled at 100 points/stride, $\tau = 10 \implies 0.8\tau = 8$ samples (8% of a single stride). Neighbor search picks points along the same instantaneous gait cycle!

**Correct Principle:**
Theiler window $W$ must satisfy:
$$W \ge T_{stride} \quad \text{or} \quad W \ge 0.5 T_{stride}$$
guaranteeing that nearest neighbors are strictly drawn from **different movement cycles** (orbits).

---

## 5. Temporal Normalization: Continuous vs. Stride-by-Stride

### Method 1: Continuous Whole-Sequence Normalization
Given $N_{strides}$ with start frame $f_1$ and end frame $f_{end}$:
$$s_{norm}(t) = \text{Interp}\left(\{s(f)\}_{f=f_1}^{f_{end}}, \; N_{samples} = N_{strides} \times 100\right)$$
- Each stride has $\approx 100$ samples on average.
- **Inter-stride temporal variability is preserved:** $\Delta t_i \ne \Delta t_j$.
- Long-term exponent $\lambda_L$ correctly reflects orbital drift ($\lambda_L > 0$).

### Method 2: Stride-by-Stride Individual Normalization
For each stride $i \in [1, N_{strides}]$:
$$s_i(t) = \text{Interp}\left(\{s(f)\}_{f=f_i}^{f_{i+1}}, \; 100 \text{ points}\right)$$
- All strides are forced to have identical duration (100 samples).
- **Inter-stride phase jitter is artificially erased.**
- The long-term divergence curve collapses, causing $\lambda_L$ to become artifactually negative ($\lambda_L < 0$).

---

## 6. Kinematic State Variables: Velocity vs. Angle

The Takens theorem assumes a stationary attractor. 
- **Joint Flexion/Extension Angles:** Exhibit postural shifts, sensor baseline wander, and non-stationarity over minutes of running. This manifests as an artificially steep initial slope ($\lambda_S \approx 2.0$).
- **Joint Angular Velocity ($\dot{\theta} = d\theta/dt$):** Inherently zero-mean, stationary, and bounded, fulfilling attractor compactness requirements.
- **Rule of Physical Differentiation:** Differentiation must be performed at raw sampling frequency $F_s$ (e.g. 200 Hz) before any time-normalization:
$$\omega(t) = \frac{\theta(t + \Delta t_{raw}) - \theta(t - \Delta t_{raw})}{2 \Delta t_{raw}}$$
