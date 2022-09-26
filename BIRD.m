function EEG_out = BIRD(EEG_in, opts)

% Function for performing blink-identified robust detrending on EEG data in EEGLAB / ERPLAB format.
% This function is designed for event-related data, and will detrend chunks of data focused on events.
% The data can then be epoched with the same function.
% Otherwise, the large epochs of detrended data will be returned.
% This function assumes that the pre-processing steps prior to epoching have been performed,
% and an EVENTLIST and BDF have already been applied to the data.
% Additionally, a VEOG channel must be created and specified.
%
% Inputs:
%    EEG: EEGLAB data structure
%    opts: A struct containing the following fields:
%        REQUIRED:
%            veog: Integer specifying the channel number of the VEOG channel
%        OPTIONAL:
%            length: Single value specifying time length of epoch for detrending. 
%                    Should be specified in seconds. (default = 35) 
%            order: Polynomial order for the robust detrending process. (default = 60)
%            threshold: Amplitude threshold for identifying blinks. (default = 120)
%            epoch: Boolean (1 or 0) specifying whether to epoch after detrending. (default = 0)
%            epochtime: 2-value vector specifying the pre-event and post-event epoch time.
%                       Same as in ERPLAB. (E.g. [-200 1000] is 200 ms pre-event, 1000 ms post-event)
%            epochbaseline: String specifying baseline option for epoching, as in ERPLAB. 
%                           Can be 'pre', 'post', 'all', or 'none'. (default = 'pre')
%            parallel: Boolean (1 or 0) specifying whether to parallelize detrending.
%                      Will loop through detrending channels using parfor.
%                      Potentially useful when number of channels is high. (default = 0)
%
% Outputs:
%    EEG: EEGLAB data structure
%
% Author: Ryan J. Hubbard, UIUC / Beckman Institute, 2022
% See Github: 
%
% BIRD uses functions from NoiseTools by Alain de Cheveigné. See:
%
% de Cheveigné, A., & Arzounian, D. (2018). 
% Robust detrending, rereferencing, outlier detection, and inpainting for multichannel data. 
% NeuroImage, 172, 903-912.
%
% BIRD uses functions from EEGLAB.
% EEGLAB is Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% BIRD uses functions from ERPLAB.
% ERPLAB is Copyright (C) 2007 The Regents of the University of California
% Created by Javier Lopez-Calderon and Steven Luck
% Center for Mind and Brain, University of California, Davis,
% javlopez@ucdavis.edu, sjluck@ucdavis.edu

if nargin < 2
    error('Two inputs - EEG and opts - required!');
end
if ~isstruct(EEG_in) || ~isstruct(opts)
    error('EEG and opts must be structures!');
end
if isempty(fieldnames(opts))
    error('opts structure is empty!');
end
if ~isfield(EEG_in, 'EVENTLIST')
    error('BIRD requires EVENTLIST to be present in EEG structure!');
end
if ~isfield(opts, 'veog')
    error('VEOG channel must be specified in opts!');
end
if ~isfield(opts, 'length') opts.length = 35; end
if ~isfield(opts, 'order') opts.order = 60; end
if ~isfield(opts, 'threshold') opts.threshold = 120; end
if ~isfield(opts, 'epoch') opts.epoch = 0; end
if ~isfield(opts, 'epochtime') opts.epochtime = [-200 1000]; end
if ~isfield(opts, 'epochbaseline') opts.epochbaseline = 'pre'; end
if ~isfield(opts, 'parallel') opts.parallel = 0; end

EEG = EEG_in;
eegbdf = EEG.EVENTLIST.bdfname;
EEG.urevent = [];
EEG2 = pop_select( EEG,'time',[0 25] );
EEG2 = pop_selectevent( EEG2, 'omitlatency','0<=25','deleteevents','on');
EEG2.data = flip(EEG2.data,2);
EEG3 = pop_mergeset( EEG2, EEG, 0);
EEG3  = pop_creabasiceventlist( EEG3 , 'AlphanumericCleaning', 'on', 'Newboundary', { -99 }, 'Stringboundary', { 'boundary' }, 'Warning', 'off' );
EEG3  = pop_binlister( EEG3 , 'BDF', eegbdf, 'ExportEL', 'no', 'ImportEL', 'no', 'Saveas', 'off', 'SendEL2', 'EEG&Text', 'Warning', 'off' );
EEG = EEG3;

EEG2 = pop_select( EEG,'time',[floor(EEG.times(end)/1000)-25 floor(EEG.times(end)/1000)] );
EEG2 = pop_selectevent( EEG2, 'omitlatency','0<=25','deleteevents','on');
EEG2.data = flip(EEG2.data,2);
EEG3 = pop_mergeset( EEG, EEG2, 0);
EEG3  = pop_creabasiceventlist( EEG3 , 'AlphanumericCleaning', 'on', 'Newboundary', { -99 }, 'Stringboundary', { 'boundary' }, 'Warning', 'off' );
EEG3  = pop_binlister( EEG3 , 'BDF', eegbdf, 'ExportEL', 'no', 'ImportEL', 'no', 'Saveas', 'off', 'SendEL2', 'EEG&Text', 'Warning', 'off' );
EEG = EEG3;
clear EEG2
clear EEG3

preepoch = EEG;

twin = (opts.length * 1000) / 2;
EEG = pop_epochbin( EEG , [(twin*-1)  twin],  'none');

dat = permute(EEG.data,[2 1 3]);
VEOG = squeeze(dat(:,opts.veog,:));

detrend_dat = zeros(size(dat));

for i = 1:size(dat,3)
    VEOG1 = VEOG(:,i)- mean(VEOG(:,i));
    VEOG1 = double(VEOG1 > opts.threshold);
    [labeledRegions, numRegions] = bwlabel(VEOG1);
    wtst = ones(size(dat,1),size(dat,2));
    if numRegions > 0
        for j = 1:numRegions
            fdfd = find(labeledRegions == j);
            tmp1 = fdfd(1) - 50;
            if tmp1 < 1
                tmp1 = 1;
            end
            tmp2 = fdfd(end) + 50;
            if tmp2 > size(wtst,1)
                tmp2 = size(wtst,1);
            end
            wtst(tmp1:tmp2,:) = 0;
        end
    end

    clear ME1

    try
        [y4,w4,r4]=nt_detrend(squeeze(dat(:,:,i)),1,wtst,'polynomials',3,10);
    catch ME1
        if strcmp(ME1.message, 'weights all zero')
            wtst = ones(size(dat,1),size(dat,2));
            [y4,w4,r4]=nt_detrend(squeeze(dat(:,:,i)),1,wtst,'polynomials',3,10);
        end
    end   

    clear ME1
    
    if ~opts.parallel

        try
            [y3,w3,r3]=nt_detrend(y4,opts.order,w4);
        catch ME1
            if strcmp(ME1.message, 'weights all zero')
                wtst = ones(size(dat,1),size(dat,2));
                [y4,w4,r4]=nt_detrend(squeeze(dat(:,:,i)),1,wtst,'polynomials',3,10);
                [y3,w3,r3]=nt_detrend(y4,opts.order,w4);
            end
        end
        
    else
        
        y3 = zeros(size(y4));
        p = gcp;
        pctRunOnAll warning('off','all')
        parfor j = 1:size(wtst,2)
            [yy,w3,r3]=nt_detrend(squeeze(y4(:,j)),60,squeeze(w4(:,j)));
            y3(:,j) = yy;
        end
        
    end

    detrend_dat(:,:,i) = y3;
end

detrend_dat = permute(detrend_dat,[2 1 3]);

if opts.epoch
    EEG = preepoch;
    EEG = pop_epochbin( EEG , opts.epochtime,  opts.epochbaseline);
    eegsr = 1000 / EEG.srate;
    detrend_dat = detrend_dat(:,((size(detrend_dat,2)/2)+1)-abs(round(opts.epochtime(1)/eegsr)):((size(detrend_dat,2)/2))+abs(round(opts.epochtime(2)/eegsr)),:);    
    EEG.data = detrend_dat;
else
    EEG.data = detrend_dat;
end

EEG_out = EEG;
end