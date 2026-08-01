# Stage A dataset validation

Stage A checks the three datasets selected by the user. It does not edit the
experiment folders, recreate triggers, filter spikes, or copy the large
`sp_corr` arrays.

## Run order

1. Copy the complete `stage_a_check` folder to any location. Keep the
   `helpers` subfolder inside it.
2. Open `00_CheckSettings.m` and set `single_folder`, `sim_folder`, and
   `seq_folder`.
3. Confirm `FS` and `Electrode_Type`. The AnalysisFunctions path is derived
   automatically from the dataset location.
4. Run `01_StageA_ValidateDatasets.m`.
5. Review the complete report printed in the MATLAB Command Window. The same
   report is also saved as `check_results/StageA_Report.txt`.

If simultaneous and sequential conditions are stored in the same experiment
folder, set `sim_folder` and `seq_folder` to the same path.

Stage A saves decoded metadata in
`check_results/StageA_DatasetInfo.mat`. Stage B will read that file after the
Stage A checks have been confirmed on real data.

## Pass conditions

For each paired dataset, Stage A expects:

- one `*sp_xia_SSD.mat` file containing `sp_corr`;
- one `*_exp_datafile_*.mat` file;
- one existing `*.trig.dat` file;
- a bad-trial file with the same global trial list repeated across channels;
- a bad-channel file; and
- a final `*_RespondingChannels.mat` file.

The single dataset is allowed to have no bad-trial, bad-channel, or responding
file. This is reported as information rather than treated as an error.
