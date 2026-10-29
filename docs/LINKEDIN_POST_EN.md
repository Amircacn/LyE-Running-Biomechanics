# Ready-to-Publish LinkedIn Post (English)

## Option A: Punchy & Concise (~1,850 characters - Recommended)
*(Fits safely within LinkedIn character limits with zero trimming)*

```markdown
Have you ever used an established algorithm from a paper or toolbox, fed it your data, gotten a number, and moved on?

In nonlinear dynamics and biomechanics, this "black-box" trust can silently invalidate your conclusions.

Recently, I conducted a computational audit of the Largest Lyapunov Exponent (LyE) and Local Dynamic Stability (LDS)—metrics widely used to quantify locomotor stability in gait and running.

I benchmarked two gold-standard reference implementations on real 200 Hz running kinematics (80 full strides):
1. The classic NONAN Toolbox code (Rosenstein et al., via Cignetti & Senderling)
2. The multi-neighbor framework by Dr. Sina Mehdizadeh (2019)

Running them head-to-head revealed three critical methodological pitfalls:

1. The Theiler Window Trap
In classic NONAN, the exclusion window was hard-coded as W = 0.8*tau. In 200 Hz running data (tau ~ 10), that window was only 8 samples—just 8% of a stride! The algorithm matched "nearest neighbors" along the very same stride, measuring velocity along its own path rather than true divergence. Decoupling W to >= 50% of a stride fixed this within-orbit bias.

2. Single-Neighbor Noise & The Stability Plateau
Classic implementations track only k = 1 neighbor, making them fragile against MoCap marker noise (lambda_S falsely jumped to 2.18). Averaging divergence across k neighbors (Mehdizadeh, 2019) filters this noise. Sweeping k = 1..30 on our knee kinematics showed divergence stabilizing into a rock-solid plateau starting at k = 15 (lambda_S = 1.57, R^2 = 0.972).

3. Center of Mass vs. Joint Angles
Most gait literature tuned pipelines on Center of Mass (CoM) velocity (lambda_S ~ 0.7-0.9). On knee joint angle, lambda_S reaches 1.98. Whole-body linear momentum acts as a low-pass filter, while segmental rotations exhibit sharp state-space curvature. Comparing them without tuning geometry is an apples-to-oranges mistake.

The Takeaway:
Nonlinear algorithms aren't plug-and-play. Always verify embedding parameters, decouple the Theiler window, and sweep k for your specific sensor.

All MATLAB code, benchmarks, and the 200 Hz dataset are open-source on GitHub:
https://github.com/Amircacn/LyE-Running-Biomechanics

How do you tackle parameter sensitivity in your biomechanical workflows?

#Biomechanics #DataScience #NonlinearDynamics #LyapunovExponent #SportsScience #MATLAB #GaitAnalysis
```
