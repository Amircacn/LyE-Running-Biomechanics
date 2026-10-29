# Ready-to-Publish LinkedIn Post (English)

> **Publication Guide:**  
> 1. Copy the text below directly into LinkedIn.
> 2. Attach `docs/figures/Master_Biomechanics_Dashboard_White.png` or `docs/figures/State_Space_Attractor_3D.png`.
> 3. The link `https://github.com/Amircacn/LyE-Running-Biomechanics` is already embedded.

---

### Why Tuning Your Computational Pipeline Matters More Than You Think in Biomechanics & Nonlinear Data Analysis

Have you ever taken an established algorithm from a published paper or toolbox, fed it your experimental data, gotten an output number, and moved on?

In nonlinear dynamics and biomechanics, this "black-box" trust can silently invalidate an entire study.

Recently, I conducted a deep computational audit of the **Largest Lyapunov Exponent (LyE)** and **Local Dynamic Stability (LDS)** algorithms—tools widely used to quantify neuromuscular stability during human locomotion (gait and running). 

I pitted two gold-standard reference implementations:
1. The classic **NONAN Toolbox** implementation (Rosenstein et al., via Cignetti & Senderling / University of Nebraska Omaha).
2. The multi-neighbor framework proposed by **Dr. Sina Mehdizadeh (2019)**.

Running them head-to-head on **real 200 Hz optical motion capture running kinematics (80 full strides)** and canonical chaotic systems revealed critical methodological pitfalls:

---

#### 1. The Theiler Window Trap in Classic Implementations
In the NONAN code, the Theiler exclusion window was hard-coded as $0.8\tau$ (a fraction of embedding delay). In 200 Hz running data (100 samples/stride, $\tau \approx 10$), that window is just **8 samples**—barely 8% of a single stride! 
The algorithm ended up picking "nearest neighbors" along the *very same stride*. Instead of measuring the divergence of neighboring trajectories in state space, it was tracking the velocity of points sliding along their own orbital trajectory! 
Decoupling the Theiler window to at least half a stride ($W \ge 50$ samples) fixed this within-orbit contamination.

#### 2. The Single-Neighbor Noise Vulnerability
Classic algorithms track only one neighbor ($k = 1$). High-frequency optical marker wobble immediately inflates short-term exponents ($\lambda_S \approx 2.18$). 
Mehdizadeh’s multi-neighbor insight averages trajectory divergence across $k$ neighbors. In our knee velocity data, sweeping $k = 1 \dots 30$ revealed an initial steep noise-filtering drop before hitting a crystal-clear, rock-solid plateau starting at **$k = 15$** ($\lambda_S = 1.5701, R^2 = 0.972$).

#### 3. Center of Mass (CoM) vs. Joint Angles: Apples and Oranges
Most foundational walking studies tuned their pipelines on **Center of Mass (CoM) linear velocity**, yielding $\lambda_S \approx 0.7 - 0.9$. 
Run the exact same tuned pipeline on **knee joint flexion angle**, and $\lambda_S$ jumps to **1.98** (or **1.57** for knee angular velocity). 
Why? Whole-body linear momentum acts as a massive low-pass filter, whereas segmental joint rotations exhibit sharp angular curvature in state space. **Comparing stability numbers across different physiological signals without tuning state-space geometry is an apples-to-oranges mistake.**

#### 4. The Time-Normalization Pitfall
Normalizing each stride individually to 100 points (stride-by-stride warping) artificially crushes inter-stride temporal variance, driving long-term divergence negative ($\lambda_L < 0$). Normalizing the entire continuous sequence preserves true neuromuscular variability ($\lambda_L = +0.058$).

---

### 💡 The Takeaway
Nonlinear time-series algorithms are not "plug-and-play." 
- Always verify your embedding parameters ($\tau$, $m$) with AMI and FNN.
- Decouple the Theiler window based on your signal's physical period.
- Sweep $k$-neighbors to locate the stability plateau for your specific sensor.
- Differentiate raw signals before resampling, never after.

To support reproducible research, I built an open-source, fully configurable, and benchmarked MATLAB framework with our 200 Hz running dataset:

🔗 **GitHub Repository:** https://github.com/Amircacn/LyE-Running-Biomechanics

I'd love to hear your thoughts: How do you handle parameter sensitivity when dealing with noisy biological time series?

#Biomechanics #DataScience #NonlinearDynamics #LyapunovExponent #SportsScience #MATLAB #GaitAnalysis #MotionAnalysis #BiomedicalEngineering #OpenScience
