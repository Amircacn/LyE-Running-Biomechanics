# Curated Biomechanical Sample Data

This directory contains verified sample data to enable instantaneous reproduction of all examples:

1. `sample_running_knee_angles.csv`:
   - Source: 200 Hz optical motion capture of treadmill running (Subject 01, Control condition).
   - Columns: `Frame`, `Time_s`, `R_Knee_FlexExt_deg`, `R_Knee_AbdAdd_deg`, `R_Knee_IntExtRot_deg`.
   - Length: 24,001 frames (120 seconds of continuous running, ~151 full strides).

2. `sample_com_velocity.mat`:
   - Source: Optical motion capture of locomotion.
   - Variables: `CoM_ML` (Medio-lateral position), `CoM_vel` (physical derivative at raw Fs), `lhs_events` (heelstrike sample indices), `fs_opto` (sampling rate = 100 Hz).
