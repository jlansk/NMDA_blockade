# Memantine and Alzheimer’s disease effects on NMDA-receptor blockade
_This code accompanies Lanskey et al. (2024)._

The code runs on Matlab 2020b and uses SPM12 v7771, which is freely available:
* [Statistic Parametric Mapping](http://www.fil.ion.ucl.ac.uk/spm/)

## Analysis steps
The analysis is run in the following order:
1. Preprocessing, using same pipeline as dcm_cmc repository, please see folder 'jlansk/dcm_cmc_ntad/scripts/1_preprocessing/'](https://github.com/jlansk/dcm_cmc_ntad/tree/main/scripts/1_preprocessing'
2. Sensor space analysis and magnesium block switch plot, please see folder in this repository: '2_sensor'
3. Part 1: First-level DCM inversion and PEB of memantine versus placebo in healthy controls, folder '3_part1'
4. Part 2: First-level DCM inversion and PEB of people with Alzheimer's disease, folder 4_part2
