# 13C-DMFA
A MATLAB-based workflow to perform dynamic flux analysis using time-course 13C stable isotope data, as well as analysis in the steady-state regime. It provides a text-based interface to control the underlying modelling architecture. The workflow is scripted and executed in stages (see Figure), but are open for moodifications (e.g., **leastSQ_.m**). Each OpenFLUX object is self-contained (data, model and parameters are uploaded), and is portable (not dependent on the original inputs).
<br />
<br />
![alt text](https://github.com/lakeeeq/OpenFLUX/blob/master/OpenFLUX%20workflow.png)
<br />
<br />
## Getting Started
### Prerequisite
MATLAB (2017 or later) with Optimisation Toolbox.

### Installation
Download all scripts/functions. Start at **OFstartHere.m**, using MATLAB comment "%" to stage workflow. Modify parameters in **OFspec_SBR.m** or **OFspec_ODE.m** or **OFspec_SS.m** to configure OpenFLUX objects.

## Running tests
1. Adipocyte model and data described in manuscript provided (folder *inputs_adipocytes*).
2. Toy dynamic model.
3. Toy steady-state model.
4. Toy steady-state model for OpenFLUX version 2009 workflow (use **OFstartHere_2009.m**).

### Main purpose of each scripts
**runner_1_modelSetup.m** to setup OpenFLUX object for first time.

**runner_2_genFeasible.m** to find a feasible solution required to simulate. No data required.

**runner_3a_optimisationSetup.m** to read data and set up OpenFLUX objects for multi-start optimisation.

**runner_3b_optimisationRun.m** or **runner_3b_optimisationRun_HPC.m** to run optimisation with data as constraints.

**runner_3c_changeStepSize.m** or **runner_3d_changeSBRtoODE.m** to modiffy OpenFLUX objects.

**runner_4a_visualiseSoln.m** to plot the feasible solution contained in an OpenFLUX object.

**runner_4b_visualiseMetData.m** to plot measured isotopologue data.

**runner_5a_MonteCarloSetup.m** to set up OpenFLUX objects with input data corrupted. Then use *runner_3b...* to optimise.

**runner_5b_compileVisualiseOps.m** to compile optimised OpenFLUX objects and visualise results.

## Contact
Lake-Ee Quek, lake-ee.quek@sydney.edu.au

## Reference
Lake-Ee Quek, James R. Krycer, Satoshi Ohno, Katsuyuki  Yugi, Daniel J. Fazakerley, Richard Scalzo, Sarah D. Elkington, Ziwei Dai, Akiyoshi Hirayama, Satsuki Ikeda, Futaba Shoji, Kumi Suzuki, Jason W. Locasale, Tomoyoshi Soga, David E. James, Shinya Kuroda. *Dynamic 13C flux analysis captures the reorganisation of adipocyte glucose metabolism in response to insulin.* Submitted.
