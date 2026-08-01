# Stage A dataset validation

Stage A checks the three datasets selected by the user. It does not edit the
experiment folders, recreate triggers, filter spikes, or copy the large
`sp_corr` arrays.

## Run order

1. Copy the complete `stage_a_check` folder to any location. Keep the
   `helpers` subfolder inside it.
2. Open `StageA_ValidateDatasets.m` and edit only the `USER SETTINGS` section:
   the three dataset paths, electrode type, and sampling rate.
3. If Sim and Seq share one folder, set `seq_folder = sim_folder`.
4. Run the complete script. It changes folders internally when needed; no
   manual `cd` and no pop-up windows are required.
5. Review the complete report printed in the MATLAB Command Window. The same
   report is also saved as `check_results/StageA_Report.txt`.

The AnalysisFunctions path is derived automatically from the Single dataset
path. If it cannot be found, the script reports the expected location. The
current package's helper folder is explicitly kept above older checking-code
copies on the MATLAB path.

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
