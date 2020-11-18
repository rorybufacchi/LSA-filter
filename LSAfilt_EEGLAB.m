function [EEG_LSA EEG_LSA_nan EEG_LSA_selection] = LSAfilt_EEGLAB(EEG, varargin)

% =========================================================================
% set default for non-specified options
if nargin == 1
    ops     = LSAfilt_setDefaults(EEG);
elseif nargin == 2
    ops     = LSAfilt_setDefaults(EEG,varargin{1});
else
    ops     = CreateOps(varargin);
    ops     = LSAfilt_setDefaults(ops);
end

% Dimensions over which to average for the wave plot
avDims = [];                
if ops.avChanFl == 1
    avDims = [avDims 1];
end
if ops.avTrialFl == 1
    avDims = [avDims 3];
end
% =========================================================================

% =========================================================================
% Perform data transformations: permute order of eeg
eegTmp = permute(EEG.data,[2 3 1]);
% =========================================================================     


% =========================================================================
% Reformat the data
pltEEG = squeeze(nanmean(EEG.data(ops.chan, ops.waveTP, ops.trl),avDims));
% Ensure the dimensions are correctly ordered
if ismember(1,avDims)
    pltEEG = permute(pltEEG,[4 1 2 3]);
end
% =========================================================================


% =========================================================================
% Define the output EEG variables
[lsaEEG lbd]            = LSAfilt(eegTmp,ops.ref,'adaptive');
filtEEGadapt            = permute(lsaEEG,[3 1 2]);
lbd                     = permute(lbd,[3 1 2]);
filtEEGstandard         = permute(LSAfilt(eegTmp,ops.ref,'standard'),[3 1 2]);

sn                      = sign(ops.scalpTPend - ops.scalpTP);
if sn == 0
    sn = 1;
end
eegSel                  = nanmean(eegTmp(ops.scalpTP:sn:ops.scalpTPend,:,:),1);
for iCh =1:EEG.nbchan
    % First dimension is signal channels, second is reference channels
    [tmpEEGadapt(:,:,:,iCh), tmp(:,iCh)]     = LSAfilt(eegSel,iCh,'adaptive');
end
filtEEGadaptTP          = squeeze(tmpEEGadapt);
filtEEGadaptTP          = permute(filtEEGadaptTP (:,:,ops.ref),[2 3 1]);
lbdTP                   = squeeze(tmp);
% =========================================================================


% =========================================================================
% Plot the interactive figure if requested
if ops.guiFl == 1
    global gEEG gOps gH gV
    
    LSAfilt_TimePlotInteractive;
    waitfor(gH.fig1);
    
    % Store important variables
    filtEEGadapt        = gV.filtEEGadapt;
    filtEEGadaptTP      = permute(gV.filtEEGadaptTP (:,:,gOps.ref),[2 3 1]);
    ops                 = gOps;
    
    % Clear all the global variables
    clear gEEG gOps gH gV
end
% =========================================================================


% =========================================================================
% Store new EEG structures - filtered without NaNs, filtered with NaNs and only selected timepoints

% Values which should be set to NaN (are not interpretable)
toNaN                   = repmat(lbd < 1,[1 1 size(EEG.data,3) size(EEG.data,4)]);
% Filter output (complete)
filtOut                 = filtEEGadapt;
filtOut(toNaN)          = NaN;

EEG_LSA                 = EEG;
EEG_LSA_nan             = EEG;
% Extract only one timepoint for the averaged data and then later replace it with the averaged data
timesMid                = ((EEG.times(ops.scalpTP) + EEG.times(ops.scalpTPend))./2)./1000;
jumpTP                  = mean(diff(EEG.times))./1000;
EEG_LSA_selection       = pop_select( EEG, 'time',[timesMid timesMid+jumpTP] );

% EEG_LSA = pop_editset(EEG, 'setname', 'only Dat 1TEST', 'run', []);
EEG_LSA.setname         = [EEG_LSA.setname ' LSA filtered'];
EEG_LSA_nan.setname     = [EEG_LSA_nan.setname ' LSA filtered - non suitable replaced with NaN'];
EEG_LSA_selection.setname = [EEG_LSA_selection.setname ' LSA filtered - selection only'];

EEG_LSA.data            = filtEEGadapt;
EEG_LSA_nan.data        = filtOut;
EEG_LSA_selection.data  = filtEEGadaptTP;
% =========================================================================


% Display warning message:
disp(' ');
LSA_EEGLAB_WarningMessageCmdWindow;
disp(' ');

end

function [ops] = CreateOps(args)
    % Creates the options structure that the LSA filtering requires

    for iN = 1:2:length(args)
        ops.(args{iN}) = args{iN+1};
    end
end
