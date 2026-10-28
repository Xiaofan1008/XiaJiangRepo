# Neurostim / Open Ephys Experiment Run Guide

This guide is for real stimulus runs on the stimulus PC with Open Ephys recording on the ephys PC. It assumes the current rig configuration is `marmolab.MU00011699` and Open Ephys is controlled through the REST API at:

```matlab
http://130.194.192.3:37497
```

## Before Starting

On the Open Ephys PC:

1. Open the Open Ephys GUI and load the correct Neuropixels acquisition graph.
2. Make sure the GUI is not already recording.
3. Manually create today's save folder:

```text
D:\data\YYYY\MM\DD
```

For example, on 2026-07-02:

```text
D:\data\2026\07\02
```

This folder must already exist. If it does not exist, Open Ephys may silently keep saving into the previous valid date folder.

On the stimulus PC / MATLAB:

1. Confirm MATLAB sees the correct date:

```matlab
datestr(now,'yyyy\\mm\\dd\\')
```

2. Confirm the rig config points Open Ephys to today's folder:

```matlab
c = marmolab.rigcfg('debug', true);
c.oephys.recordDir
```

Expected output for 2026-07-02:

```matlab
'D:\data\2026\07\02\'
```

3. Confirm Open Ephys API communication:

```matlab
webread('http://130.194.192.3:37497/api/status')
```

Expected output is a struct with `mode` equal to `IDLE`, `ACQUIRE`, or similar. If this errors, MATLAB cannot talk to Open Ephys.

4. Confirm key MATLAB code is on the path:

```matlab
which freeviewing.runTranslatingImages -all
which mapping.noisegridPassive -all
which marmolab.rigcfg -all
which neurostim.cic -all
which Screen -all
```

## Important Rules

- Do not use `PATH_TEST` for real recordings.
- Do not manually press Record in Open Ephys for these runs. The Neurostim/Open Ephys plugin starts recording automatically when the stimulus run begins.
- Use `debug`, `false` for real recordings.
- Create a fresh run from the function/script. Do not reuse an old `c` object from a previous day or previous session.
- The Open Ephys GUI parent-directory display may still show an old value before the run. The important check is the API state after MATLAB configures Open Ephys.

## Running Translating Images

The underlying function is:

```matlab
freeviewing.runTranslatingImages(filenames, ...
    'subject', 'SUBJECT_ID', ...
    'paradigm', 'PARADIGM_NAME', ...
    'moveSpeeds', 3, ...
    'contrasts', {1}, ...
    'debug', false);
```

`runTranslatingImages` now accepts an optional `paradigm` argument. If omitted, it defaults to `translatingImages`, but real experiment scripts should pass an explicit section-specific paradigm name.

The Neurostim file will therefore be named like:

```text
SUBJECT_ID.PARADIGM_NAME.HHMMSS.mat
```

The Open Ephys folder will be named with underscores because the Open Ephys REST API does not accept periods in `base_text`:

```text
SUBJECT_ID_PARADIGM_NAME_HHMMSS_YYYY-MM-DD_HH-MM-SS
```

## Running Sections From David_stimuli.m

Open:

```text
/home/stimadmin/Documents/MATLAB/marmolab/marmolab-stimuli/+freeviewing/David_stimuli.m
```

Run one section at a time.

Before running sections, set the real recording configuration at the top of the file:

```matlab
SUBJECT_ID = 'CJ270';

MOVE_SPEED = 3;
CONTRASTS = {1};
```

`SUBJECT_ID` controls the subject portion of both the Neurostim `.mat` filename and the Open Ephys folder name. Use a real subject ID for real recordings, and a clearly test-only ID for dry runs.

`MOVE_SPEED` controls image motion for all sections:

```matlab
MOVE_SPEED = 0;  % stationary images
MOVE_SPEED = 3;  % translating/jittering images
```

`CONTRASTS` controls image contrast over each trial:

```matlab
CONTRASTS = {1};            % constant full contrast
CONTRASTS = {[1 0 1 0]};    % flash on/off
```

Each section passes `SUBJECT_ID`, `MOVE_SPEED`, `CONTRASTS`, and its own paradigm name into `freeviewing.runTranslatingImages(...)`. The current section paradigms are:

```text
translatingImages_origRecon
translatingImages_origPhaseScrambled
translatingImages_latentWalk7
translatingImages_latentWalk13
translatingImages_latentWalk25
translatingImages_imageToImageWalk6
translatingImages_imageToImageWalk12
translatingImages_imageToImageWalk24
translatingImages_origPSVariants
translatingImages_psVariantsScrambled
```

For example, with `SUBJECT_ID = 'CJ270'`, the first section saves:

```text
Neurostim:
CJ270.translatingImages_origRecon.HHMMSS.mat

Open Ephys:
CJ270_translatingImages_origRecon_HHMMSS_YYYY-MM-DD_HH-MM-SS
```

For dry-run section testing, temporarily reduce the section's image count or folder count, then restore it before real recording.

## Running Noise Grid Passive

Use:

```matlab
mapping.noisegridPassive('SUBJECT_ID', ...
    'debug', false);
```

For example:

```matlab
mapping.noisegridPassive('CJ270', ...
    'debug', false);
```

`noisegridPassive` sets:

```matlab
c.paradigm = 'noisegrid';
```

The Neurostim file will therefore be named like:

```text
SUBJECT_ID.noisegrid.HHMMSS.mat
```

The Open Ephys folder will be named like:

```text
SUBJECT_ID_noisegrid_HHMMSS_YYYY-MM-DD_HH-MM-SS
```

## During The Run

When the stimulus run starts, Neurostim should:

1. Create the local Neurostim `.mat` file under:

```text
/home/stimadmin/data/YYYY/MM/DD/
```

2. Send Open Ephys the recording directory and session name.
3. Put Open Ephys into `RECORD`.
4. Display the stimuli on the external monitor.

Avoid typing in the MATLAB command window while the experiment is running unless you intend to interrupt/abort the run.

## After The Run

Check the Open Ephys API state:

```matlab
r = webread('http://130.194.192.3:37497/api/recording');
r.parent_directory
r.record_nodes.parent_directory
r.base_text
```

For a run on 2026-07-02, the parent directories should be:

```text
D:\data\2026\07\02
```

Check the Open Ephys PC for a new folder under:

```text
D:\data\YYYY\MM\DD\
```

Check the stimulus PC for the matching Neurostim file under:

```text
/home/stimadmin/data/YYYY/MM/DD/
```

The names will not be character-for-character identical because Neurostim uses periods and Open Ephys uses underscores. They should still match by subject, paradigm, and Neurostim `HHMMSS` session time.

Example pair:

```text
Neurostim:
CJ270.translatingImages_origRecon.143509.mat

Open Ephys:
CJ270_translatingImages_origRecon_143509_2026-07-02_14-35-09
```

## Troubleshooting Quick Checks

If Open Ephys saves to yesterday's folder:

- Manually create today's folder on the Open Ephys PC.
- Recreate the MATLAB `cic` object by starting a new run.
- Confirm `c.oephys.recordDir` points to today's date.

If Open Ephys creates `experiment2`, `experiment3`, etc. inside an old folder:

- Stop recording.
- Make sure Open Ephys is idle.
- Restart Open Ephys if necessary.
- Confirm the folder name is being accepted by the API:

```matlab
r = webread('http://130.194.192.3:37497/api/recording');
r.base_text
```

If MATLAB cannot find the run functions after restart:

```matlab
addpath(genpath('/home/stimadmin/Documents/neurostim'));
addpath(genpath('/home/stimadmin/Documents/MATLAB/marmolab/marmolab-common'));
addpath(genpath('/home/stimadmin/Documents/MATLAB/marmolab/marmolab-stimuli'));
addpath(genpath('/home/stimadmin/Documents/MATLAB/marmolab/marmodata'));
addpath('/home/stimadmin/Documents/MATLAB-git');
addpath(genpath('/home/stimadmin/Documents/Psychtoolbox-3/Psychtoolbox'));
addpath(genpath('/home/stimadmin/Documents/open-ephys-matlab-tools'));
addpath(genpath('/home/stimadmin/Documents/V1_1_urlread2'));
```

Then rerun the `which ... -all` checks above.
