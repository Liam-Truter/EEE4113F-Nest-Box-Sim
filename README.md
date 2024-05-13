# Nest Box Simulator
## Needed Packages
- numpy
- scipy
- matplotlib

## Breakdown of tests:
### Panel Test A:
Test of thermal resistances of 2cm wood panels from various databases.
### Panel Test B:
Test of thermal resistances of added insulation layers.
### Nest Box Test A:
Test of thermal performance of various panel configurations on an average summer day.
### Nest Box Test B:
Test of thermal performance of various panel configurations on an artificially hot day.
### Nest Box Test C:
Test of thermal and energy performance of various Peltier intensities on an artificially hot day.

## Known Issues
### Temperature oscillations/checkerboarding
Can occur when using a small mesh size and a larger time step. Can be fixed by increasing time resolution.
