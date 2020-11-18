%  pop_LSAfilt_EEGLAB() - Apply LSA Filter
%
% Usage:
%   EEG_LSA, com, ops, EEG_LSA_nan, EEG_LSA_selection] = pop_LSAfilt_EEGLAB(EEG); % launch a GUI
%   EEG_LSA, com, ops, EEG_LSA_nan, EEG_LSA_selection]= pop_LSAfilt_EEGLAB(EEG,'key', val);
%   EEG_LSA, com, ops, EEG_LSA_nan, EEG_LSA_selection] = pop_LSAfilt_EEGLAB(EEG, options); % A structure including all options
%
% Inputs:
%   EEG                     - Input dataset
%
% Optional Inputs
%   'ref'                   - Reference electrode number
%   'scalpTP'               - Reference timepoint
%
% Outputs:
%   EEG_LSA: This dataset contains the full filtered timecourse. 
%        IMPORTANTLY, there is no guarantee that each timepoint is interpretable. 
%
%   com: Command for EEGLAB history
%
%   ops: Input options for the LSAfilt_EEGLAB function can also be input
%        into pop_LSAfilt_EEGLAB instead of key and val
%
%   EEG_LSA_nan: This dataset contains only the interpretable parts of the filtered timecourse. 
%       All timepoints for which the assumptions do not hold are set to NaN. 
%       Note that you might not be able to plot this output in EEGLAB due to its handling of NaNs. 
%
%   EEG_LSA_selection: This dataset contains the average of the data within the time window selected by the user. 
%       Hence, this dataset contains only 1 timepoint. All channels for which the assumptions do not hold are set to NaN. Therefore, the 
%       output data are fully interpretable. Note that plotting this output in EEGLAB becomes difficult as there is only 1 timepoint.
%
% See also:
%   POP_LSAFILT_EEGLAB, LSAFILT, LSAFILT_TIMEPLOTINTERACTIVE, EEGLAB
%
% Authors: Rory Bufacchi, Cesare Magri

function [EEG_LSA, com, ops, EEG_LSA_nan, EEG_LSA_selection] = pop_LSAfilt_EEGLAB(EEG, varargin)

if nargin < 1
    help pop_LSAfilt_EEGLAB;
    return;
end
    
if isempty(EEG(1).data)
    disp('Pop_LSAfilt_EEGLAB error: cannot process empty dataset'); return;
end

% =========================================================================
% set default options
if exist('ops','var')
    [ops] = LSAfilt_setDefaults(EEG,ops);
else
    [ops] = LSAfilt_setDefaults(EEG);
end
% =========================================================================

if nargin < 2
    
    % =========================================================================
    % Make the input GUI
    res = inputgui( 'geometry', { [.5 .5 .25] [.5 .5 .25] [.5 .5 .25]}, ...
        'geomvert', [1.5 2 1.5], 'uilist', { ...
        { 'Style', 'text', 'string', 'Run the LSA GUI? (Y/N)', 'fontweight', 'bold' } ...
        { 'Style', 'popupmenu', 'string', 'Yes|No' 'tag' 'GUI'},{} ...
        ...
        { 'Style', 'text', 'string', 'Reference Channel', 'fontweight', 'bold'  } ...
        { 'Style', 'edit', 'string', EEG.chanlocs(ops.ref).labels 'tag' 'chan'}, ...
        { 'style' 'pushbutton' 'string'  '...', 'enable' fastif(isempty(EEG(1).chanlocs), 'off', 'on') ...
        'callback' 'tmpchanlocs = EEG(1).chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels}, ''withindex'', ''on'', ''selectionmode'', ''single''); set(findobj(gcbf, ''tag'', ''chan''), ''string'',tmpval); clear tmp tmpchanlocs tmpval' }, ...
        ...
        { 'Style', 'text', 'string', 'LSA scalp map timepoints (ms) [start stop]', 'fontweight', 'bold' } ...
        { 'Style', 'edit', 'string', num2str(EEG.times([floor(end/2) floor(end/2 + end/4)])) 'tag' 'scalpTimes'},{} ...
        } );
    % =========================================================================
    
    % =========================================================================
    % Store the inputs
    ops.guiFl   = mod(res{1},2) ;
    allChLab    = {EEG.chanlocs.labels};
    ops.ref     = find(strcmp(res{2},allChLab));
    
    scalpTs             = eval([ '[' res{3} ']' ]);
    [~, ops.scalpTP]    = min(abs(EEG.times-scalpTs(1)));
    [~, ops.scalpTPend] = min(abs(EEG.times-scalpTs(2)));
    avDims = [];                % Dimensions over which to average for the wave plot
    if      ops.avChanFl == 1
        avDims  = [avDims 1];
    end; if ops.avTrialFl == 1
        avDims  = [avDims 3];
    end
    % =========================================================================
    
% If options are specified, use them
elseif nargin == 2
    ops     = LSAfilt_setDefaults(EEG,varargin{1});
    
% If arguments are specified, convert them to options
else
    ops     = CreateOps(varargin);
    ops     = LSAfilt_setDefaults(ops);
end

% =========================================================================
% generate command
if nargout > 1
    args = CreateArgs(ops);
    com = sprintf('EEG = pop_LSAfilt_EEGLAB( EEG, %s);', vararg2str(args));
end
% =========================================================================


% =========================================================================
% Possibly open GUI and run the filter
LSA_EEGLAB_WarningMessagePopup;
[EEG_LSA EEG_LSA_nan EEG_LSA_selection] = LSAfilt_EEGLAB(EEG, ops);
% =========================================================================

end


function [args] = CreateArgs(ops)
    % Creates the arguments cell array that the pop_functions require

    args = {};
    allN = fields(ops);
    allR = struct2cell(ops);
    for iN = 1:length(allN)
        args = { args{:} , allN{iN}, allR{iN}};
    end
end

function [ops] = CreateOps(args)
    % Creates the options structure that the LSA filtering requires

    for iN = 1:2:length(args)
        ops.(args{iN}) = args{iN+1};
    end
end
