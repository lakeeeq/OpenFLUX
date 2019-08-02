# 13C-DMFA
A MATLAB-based workflow to perform dynamic flux analysis using time-course 13C stable isotope data. It provides a text-based interface to control the underlying modelling architecture. The workflow is scripted and executed in stages (see Figure), but are open for moodifications. Each OpenFLUX object is self-contained (data, model and parameters are uploaded), and is portable (not dependent on the original inputs).

## Getting Started
### Prerequisite
MATLAB (2017 or later) with Optimisation Toolbox.

### Installation
Download all scripts/functions. Start at *OFstartHere.m*, using MATLAB comment "%" to stage workflow. Modify parameters in *OFspec_SBR.m* or *OFspec_ODE.m* to customise work.

## Running tests
Adipocyte model and data described in manuscript provided in the tests. Toy ODE model (to be included).

### Main purpose of each scripts
*runner_1_modelSetup.m* to setup OpenFLUX object for first time.
*runner_2_genFeasible.m* to find a feasible solution required to simulate.
*runner_3a_optimisationSetup.m* to set up OpenFLUX objects for multi-start optimisation.
*runner_3b_optimisationRun.m* or *runner_3b_optimisationRun_HPC.m* to run optimisation.
*runner_3c_changeStepSize.m* or *runner_3d_changeSBRtoODE.m* to modiffy OpenFLUX objects.
*runner_4a_visualiseSoln.m* to plot feasible solution contained in an OpenFLUX object
*runner_4b_visualiseMetData*.m to plot measured isotopologue data
*runner_5_MonteCarloSetup.m* to set up OpenFLUX objects with input data corrupted

## Contact
Lake-Ee Quek, lake-ee.quek@sydney.edu.au

## Reference
Lake-Ee Quek, James R. Krycer, Satoshi Ohno, Katsuyuki  Yugi, Daniel J. Fazakerley, Richard Scalzo, Sarah D. Elkington, Ziwei Dai, Akiyoshi Hirayama, Satsuki Ikeda, Futaba Shoji, Kumi Suzuki, Jason W. Locasale, Tomoyoshi Soga, David E. James, Shinya Kuroda. *Dynamic 13C flux analysis captures the reorganisation of adipocyte glucose metabolism in response to insulin.* Submitted.
