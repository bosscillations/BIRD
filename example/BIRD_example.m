% Example code for running Blink-Identified Robust Detrending on EEG data
% Note that it is suggested to low-pass filter your data before running BIRD - the example here shows the minimal steps required

% First, open EEGLAB and load in example dataset. This assumes the dataset and bdf are in the current MATLAB directory
eeglab;
currentfolder = pwd;
EEG = pop_loadset('filename', ['BIRD_data.set'], 'filepath', currentfolder);

% The example dataset sampling rate is 1000 Hz. We downsample to 250 Hz to make things faster
EEG = pop_resample(EEG,250);

% We now create the VEOG channel, which is critical for blink detection
EEG = pop_eegchanoperator( EEG, {'ch31 = ch3 - ch5 Label VEOG'});

% We now create the EVENTLIST and assign bins
EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'Eventlist', ['event.txt'], 'Newboundary', { -99 }, 'Stringboundary', { 'boundary' }, 'Warning', 'off' );
EEG  = pop_binlister( EEG , 'BDF', 'BIRD_bdf.txt', 'ExportEL', ['event_binned.txt'], 'ImportEL', 'no', 'Saveas', 'off', 'SendEL2', 'EEG&Text', 'Warning', 'off' );

% We now create the opts structure for running BIRD
opts = struct;
opts.veog = 31; % Specifies the VEOG channel number
opts.length = 30; % Specifies the length of epoch time for detrending
opts.order = 60; % Specifies the polynomial order for detrending
opts.threshold = 120; % Specifies the amplitude threshold for blink detection
opts.epoch = 1; % Specifies whether to epoch after detrending or not
opts.epochtime = [-200 1000]; % If epoching, specifies the epoch time window
opts.epochbaseline = 'pre'; % If epoching, specifies the baseline method
opts.parallel = 0; % Specifies whether to run in parallel over channels

% Finally, we run BIRD on our data!
EEG = BIRD(EEG,opts);