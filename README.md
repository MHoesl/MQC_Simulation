# MQC Simulation

Multi quantum coherences (MQC) simulation, simulates Phase Cycle Choices for a three-pulses sequence, typically used in Sodium TQ acquisitions in Nucelar Magnetic Rensonance (NMR) and Magnetic Resonance Imaging (MRI).


Tested with Matlab_R2017


## Explore the different physe cycle options with the Matlab scripts:

- ExploreTQTPPICycle.m, ExploreTQTPPICycleRelaxation.m
- ExploreFleysherRelaxation.m
- ExploreSistinaCycle.m

The scripts are using the functions:
- Coherence_Nasignal.m
- Coherence_NasignalRelaxation.m

## Results:

Phase Cycle Signal along echo Time and along B0 inhomogeneity can be simulated and the Spectra are obtained after Fourier Transform.
An exemplary result can be seen in the figure ExampleResult_FleysherSimulation.png

For visualization the following toolbox is nice and the function as.m was used: https://github.com/tsumpf/arrShow

