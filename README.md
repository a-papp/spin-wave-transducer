# Spin-wave transducer circuit model calculation
We created a circuit model of spin-wave transducers based on impedance matrices. 
The model contains three main steps:
1) Electromagnetic simulation of the current distribution and magnetic field in the tranducer lines. This is done by FEMM (https://www.femm.info/wiki/HomePage).
2) Micromagnetic simulation of the spin-waves excited by the transducers' magnetic fields. This is done by mumax3 (https://mumax.github.io/).
3) The impedance parameters (Z matrix) are calculated using a custom MATLAB script from the magnetization dynamics.

For more information, check out our paper: https://arxiv.org/abs/2410.14370
