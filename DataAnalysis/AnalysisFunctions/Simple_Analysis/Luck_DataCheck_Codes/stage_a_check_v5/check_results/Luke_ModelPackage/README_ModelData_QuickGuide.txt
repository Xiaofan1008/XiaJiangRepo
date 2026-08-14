GUIDE TO THE STIMULATION AND SPIKE DATA
=======================================

The MATLAB file contains one main structure called ModelData.

Load it using:

    S = load('DX014_D4_Pair_A030_A031_ModelData.mat');
    D = S.ModelData;


1. STIMULATION CONDITIONS
-------------------------

The five stimulation conditions are stored in:

    D.conditions

They always have the following indices:

    Index  Code       Stimulation
    -----  ---------  ------------------------------------------
      1    A          Electrode A alone at 0 ms
      2    B          Electrode B alone at 0 ms
      3    AB         A and B simultaneously at 0 ms
      4    A_to_B     A at 0 ms, followed by B at 5 ms
      5    B_to_A     B at 0 ms, followed by A at 5 ms

The stimulation electrodes are stored in:

    D.stimulation.electrode_A
    D.stimulation.electrode_B

The distance between them is stored in:

    D.stimulation.pair_distance_um

To display all condition names:

    {D.conditions.code}
    {D.conditions.label}

For a selected condition, its electrode order and pulse times are:

    D.conditions(condition_index).electrode_order
    D.conditions(condition_index).pulse_times_ms


2. STIMULATION AMPLITUDES
-------------------------

The available amplitudes are stored in:

    D.amplitudes_uA

For this dataset, they are:

    [5 6 10]    % microamps

Find the index for a particular amplitude using:

    amplitude_uA = 5;
    amplitude_index = find(abs(D.amplitudes_uA-amplitude_uA) < 1e-6,1);

The selected condition and amplitude can then be opened as one block:

    block = D.conditions(condition_index).amplitude(amplitude_index);

The amplitude is also recorded inside the block:

    block.amplitude_uA


3. SPIKE TIMES
--------------

The spike times are stored at:

    D.conditions(condition_index) ...
     .amplitude(amplitude_index) ...
     .spike_times_ms{trial_index,channel_index}

The data organization is:

    condition
      -> amplitude
        -> trial
          -> recording channel
            -> spike times

spike_times_ms is a cell array with:

    rows    = trials
    columns = recording channels

Each cell contains all spike times for one trial and one recording channel.

Example: obtain the spikes from A+B simultaneous stimulation, 5 microamps,
trial 3, and package channel 7:

    condition_index = 3;  % A+B simultaneous
    amplitude_index = find(abs(D.amplitudes_uA-5) < 1e-6,1);

    block = D.conditions(condition_index).amplitude(amplitude_index);
    spike_times = block.spike_times_ms{3,7};

All spike times are in milliseconds relative to the first stimulation pulse:

    0 ms          = first stimulation pulse
    5 ms          = second pulse for sequential stimulation
    negative time = before stimulation
    positive time = after stimulation

The exported spike-time range is stored in:

    D.metadata.stored_spike_window_ms

For this dataset, it is:

    [-50 80] ms


4. RECORDING CHANNELS
---------------------

Information about the exported recording channels is stored in:

    D.channels

The channel_index used in spike_times_ms corresponds to:

    D.channels.ChannelIndex

The associated anatomical/depth channel is:

    D.channels.DepthChannel

For example, inspect package channel 7:

    D.channels(7,:)

To find the package channel index for depth channel 7:

    channel_index = find(D.channels.DepthChannel == 7,1);

The table also contains the distance from each recording channel to the
two stimulation electrodes:

    D.channels.DistanceToA_um
    D.channels.DistanceToB_um
    D.channels.MinimumStimDistance_um


5. TRIAL INFORMATION
--------------------

The original trial IDs for a condition and amplitude are stored in:

    block.source_trial_ids

The original source trial corresponding to row 3 of spike_times_ms is:

    source_trial_id = block.source_trial_ids(3);

Use every available trial with:

    trial_indices = block.all_trial_indices;

Use the optional equal-trial-number subset with:

    trial_indices = block.balanced_trial_indices;

Extract the selected trials and all recording channels with:

    selected_spikes = block.spike_times_ms(trial_indices,:);

The balanced subset does not remove data. All included trials remain in
block.spike_times_ms.


6. DETAILED STIMULATION PARAMETERS
----------------------------------

The names of the recorded stimulation parameters are stored in:

    D.stimulation.parameter_names

The parameters for a selected condition and amplitude are stored in:

    block.stimulation_events

For example, obtain the information for the first stimulation pulse:

    event = block.stimulation_events(1);

    event.electrode
    event.pulse_time_ms
    event.representative_parameter_row

The entries in representative_parameter_row correspond, in order, to:

    D.stimulation.parameter_names

For a two-pulse condition, the second pulse parameters are stored in:

    block.stimulation_events(2)

These fields include pulse polarity, phase duration, amplitude, pulse shape,
amplifier settling, and charge-recovery settings.


7. MINIMAL WORKING EXAMPLE
--------------------------

The following example extracts all balanced-trial spikes from depth channel
7 during A-to-B stimulation at 5 microamps:

    S = load('DX014_D4_Pair_A030_A031_ModelData.mat');
    D = S.ModelData;

    % Select condition
    condition_index = find(strcmp({D.conditions.code},'A_to_B'),1);

    % Select amplitude
    amplitude_index = find(abs(D.amplitudes_uA-5) < 1e-6,1);

    % Select recording channel
    channel_index = find(D.channels.DepthChannel == 7,1);

    % Open the selected condition/amplitude block
    block = D.conditions(condition_index).amplitude(amplitude_index);

    % Use the balanced trial subset
    trial_indices = block.balanced_trial_indices;

    % Extract spikes: one cell for each selected trial
    spike_times_by_trial = ...
        block.spike_times_ms(trial_indices,channel_index);


8. SUPPLIED PLOTTING SCRIPTS
----------------------------

The following standalone scripts can be used to visualize the data:

    StageD2_PlotRasterPSTH.m
    StageD3_PlotSpikeCountCurves.m

