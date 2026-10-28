# Add Configurable David Stimuli Paradigms And Display Settings

## Summary

This update makes `David_stimuli.m` safer and clearer for running section-by-section translating-image experiments.

It adds explicit subject/paradigm naming to each section, shares display settings across all sections, and updates `runTranslatingImages.m` so wrapper scripts can pass section-specific paradigm names instead of every run saving as `translatingImages`.

## Main Changes

- Add an optional `paradigm` argument to `freeviewing.runTranslatingImages`.
- Set `c.paradigm = args.paradigm` instead of hard-coding `translatingImages`.
- Add a guard to prevent real recordings from silently saving under the default `demo` subject.
- Add top-level `David_stimuli.m` configuration:

  ```matlab
  SUBJECT_ID = '...';
  MOVE_SPEED = 3;
  CONTRASTS = {1};
  ```

- Update each `David_stimuli.m` section to pass:

  ```matlab
  'subject', SUBJECT_ID
  'paradigm', '<section-specific paradigm>'
  'moveSpeeds', MOVE_SPEED
  'contrasts', CONTRASTS
  'debug', false
  ```

- Correct the `originals vs recons` filename logic using the proven version from `defineTranslatingImages.m`.
- Update the `originals vs PS variants SCRAMBLED` filename construction to match the previously used `defineTranslatingImages.m` logic.

## Section Paradigm Names

The current section-specific paradigms are:

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

## Expected Naming

Example Neurostim file:

```text
CJ270.translatingImages_origRecon.143509.mat
```

Example Open Ephys folder:

```text
CJ270_translatingImages_origRecon_143509_2026-07-02_14-35-09
```

## Validation Performed

- Ran section-by-section dry-run tests using real images.
- Confirmed stimuli display on the VIEWPixx stimulus monitor under X11.
- Confirmed Open Ephys recording starts/stops automatically.
- Confirmed section-specific paradigm names propagate to both Neurostim `.mat` filenames and Open Ephys folder names.

