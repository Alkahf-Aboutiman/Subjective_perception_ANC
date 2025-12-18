Active Noise Control (ANC) Subjective Evaluation – R Code

This repository provides the R code used to analyze subjective evaluations of Active Noise Control (ANC) performance in a simulated vehicle interior. The study investigates how individuals perceive the effectiveness of different ANC algorithms under varying noise conditions.

Study Design:
Participants were exposed to 27 stimuli combining:

Noise types: Motorcycle, Street, Train

Sound pressure levels: 55, 65, 72 dB(A)

ANC conditions: No control, NLMS (Normalized Least-Mean-Square), SFANC-NLMS (Hybrid Selective Fixed-Filter ANC normalized LMS)

Metrics Collected:

Perceived Annoyance (PAY)

Perceived Affective Quality (PAQ)

Perceived Loudness (PLN)

Features of the Repository:

Data preprocessing and cleaning

Non-parametric repeated-measures PERMANOVA (3-way)

Post-hoc pairwise comparisons with FDR correction

Visualization of Pleasantness vs. Eventfulness maps

Code structured for reproducibility and adaptation to similar ANC evaluation studies

Citation:
If this code is useful for your research, please cite:
@article{aboutiman2025subjective,
  title={Subjective perception analysis of active noise control algorithms in an encapsulated structure: An experimental study},
  author={Aboutiman, Alkahf and Rachman, Zulfi and Oberman, Tin and Alletta, Francesco and Kang, Jian and Karimi, Hamid Reza and Ripamonti, Francesco},
  journal={Applied Acoustics},
  volume={239},
  pages={110823},
  year={2025},
  publisher={Elsevier}
}

This repository enables researchers to replicate the analysis, explore the impact of noise type, noise level, and ANC algorithms on subjective perception, and visualize affective responses in vehicle noise environments.
