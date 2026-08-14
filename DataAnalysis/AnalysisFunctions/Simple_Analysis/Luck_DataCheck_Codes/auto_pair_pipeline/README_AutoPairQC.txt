AUTOMATIC PAIR QC PIPELINE
==========================

Run only:

    Run_AutoPairQC.m

Keep the following items together in the same folder:

    Run_AutoPairQC.m
    Pipeline_StageA.m
    Pipeline_StageB.m
    Pipeline_StageC1.m
    Pipeline_StageC2.m
    Pipeline_StageD1.m
    helpers/

Do not run the Pipeline_Stage files directly.


SETTINGS TO CHANGE
------------------

At the beginning of Run_AutoPairQC.m, enter:

    single_folder
    sim_folder
    seq_folder
    Electrode_Type
    FS

If simultaneous and sequential trials are in the same folder, set:

    seq_folder = sim_folder;

The stimulation electrodes do not need to be entered. The code discovers
all candidate pairs and processes every eligible responsive pair.


DEFAULT OUTPUT
--------------

The default result folder is:

    auto_pair_results/

The pipeline saves:

    Stage A validation result and report
    Stage B automatic pair summary and report
    One Stage C1 QC data file and report for each complete pair
    One self-contained ModelData MAT file for each complete pair
    One final AutoPairQC run-summary MAT file

ModelData files are saved under:

    auto_pair_results/Luke_ModelPackage/

These intermediate files provide an audit trail, but they do not need to be
selected manually. Original experiment files are never modified.


DEFAULT FIGURES
---------------

The default setting is:

    Plot_Figures = [1 2 5];

This displays:

    1 = trial-count overview
    2 = raster and PSTH for each selected recording channel
    5 = average amplitude-response curves

The more numerous per-trial dot and individual amplitude-curve figures can
be enabled using:

    Plot_Figures = [1 2 3 4 5];

To run validation and prepare pair QC data without displaying figures:

    Plot_QC_Figures = false;

Model-data export is enabled by default:

    Export_ModelData = true;

Existing ModelData files are not overwritten unless you deliberately set:

    allow_modeldata_overwrite = true;


AUTOMATIC PAIR SELECTION
------------------------

A pair is processed only when the paired-data response union is nonempty
and the following core conditions are structurally available at the same
amplitudes:

    A alone
    B alone
    A+B simultaneous at 0 ms
    at least one sequential order at 5 ms

Therefore, all three of these layouts are accepted:

    A, B, A+B, A->B, B->A
    A, B, A+B, A->B
    A, B, A+B, B->A

If only one sequential order is present, C1 figures and the exported
ModelData file contain four conditions. No empty reverse-order condition is
invented. If the available order changes between amplitudes and there is no
one sequential order shared by every amplitude, the pair is reported for
review instead of being silently combined.

The pipeline reports incomplete or problematic pairs but does not process
them. It never selects a "best" pair based on response size.
