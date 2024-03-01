# Memantine and Alzheimer’s disease effects on NMDA-receptor blockade
_This code accompanies Lanskey et al. (2024)._

The code runs on Matlab 2020b and uses SPM12 v7771, which is freely available:
* [Statistic Parametric Mapping](http://www.fil.ion.ucl.ac.uk/spm/)

## Analysis steps
The analysis is run in the following order:
1. Preprocessing, using same pipeline as dcm_cmc, see folder 'dcm_cmc_ntad/scripts/1_preprocessing/1_preprocessing'
2. Sensor space analysis, folder '2_sensor_space'
3. Part 1: First-level DCM inversion and PEB of memantine versus placebo in healthy controls, folder '3_part1'
4. Second-level Parameteric Empirical Bayes, folder 4_part2
