---
layout: post
title: Superconducting Dark matter detectors
---
Here is my attempt at explaining the SuperCDMS detection principle without too much technical jargon.

The use of superconducting detectors for dark matter searches is motivated by their high sensitivity. Semiconductors have a small bandgap (Si: 1.1 eV, Ge: 0.74 eV) which means that the amount of energy required to excite an electron from the valence band to the conduction band of the semiconductor is small. The detectors are often cooled and held at about 10-20 mK temperature where the sensors which collect energy are superconducting. Such extremely low temperatures are achieved with the help of dilution fridges with successive stages for cooling and circulation of a mix of He3/He4 to maintain the cold temperatures.

![Dilution fridge]({{ 'assets/img/Fridge_inside.png' | relative_url}})

The basic principle of all SuperCDMS detectors is simple: measure low energy scatters of dark matter particles in semiconductor crystals. Dark matter particles (and unfortunately, neutrons) scatter off nuclei causing Nuclear Recoils (NRs). Most other background radiation scatters off electrons causing Electron Recoils (ERs). Both types of recoil have two consequences in cryogenically cooled crystals: (1) production of minute lattice vibrations, i.e. phonons, and (2) ionization (refer figure below).

![Detection mechanism cartoon]({{ 'assets/img/Detector_concept.png' | relative_url}})

SuperCDMS detectors are equipped with tungsten Transition Edge Sensors (TESs). In superconductors, there is a critical temperature below which the material becomes superconducting. Tungsten is specfically chosen because of its tunable critical temperature to improve sensor sensitivity. Superconductors have a transition curve which shows the change in electrical resistance of the material as a function of temperature:

![Example TES transition curve]({{ 'assets/img/TES_transition_curve.png' | relative_url}})

The TES is first cooled to its superconducting state, then heated via an external electrical current to hold it about half-way on the slope between superconducting and normal state. So a partical interaction depositing some energy whill heat the TES further to push it up the transition slope. So temperature goes up -> resistance goes up -> current through the TES drops. This drop in current is what is measured as a signal. The phonon signal collected this way goes through an amplifier circuit containing [Superconducting QUantum Interference Devices (SQUIDs)](https://en.wikipedia.org/wiki/SQUID) before reaching the trigger system. When the amount of energy collected breaches a certain preset threshold, the hardware issues a trigger and makes a record of the waveform.