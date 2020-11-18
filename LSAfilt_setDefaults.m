function [opsNew] = LSAfilt_setDefaults(EEG,varargin)
%[opsNew] = LSAfilt_setDefaults(ops)
%   set the default options for those which haven't been specified by ops

% Filtering options
ops.ref         = 17;                   % Reference electrode
ops.mtd         = 'adaptive';           % Filtering method

% plotting options
ops.chan        = 1:size(EEG.data,1);   % Channels to plot
ops.scalpTP     = 1;                    % Timepoints at which to plot scalp maps
ops.scalpTPend  = 1;                    % Timepoints at which to plot scalp maps
ops.waveTP      = 1:size(EEG.data,2);   % Timepoints at which to timecourse
ops.trl         = 1:size(EEG.data,3);   % Trials to plot
ops.avChanFl    = 0;                    % Flag to set averaging of channels in plot
ops.avTrialFl   = 1;                    % Flag to set averaging of trials in plot
ops.guiFl       = 1;                    % Flag to run the GUI

% Heuristic options
ops.lbdThr      = [0.2 1.0];              % Lambda limits for which to show LSA is applicable

if numel(varargin)>0
    opsNew = varargin{1};
    
    opsFields = fieldnames(ops);
    
    for iFN = 1:length(opsFields)  %field name counter
        
        % If the substructure doesn't exist, put it there in its entirety
        if ~isfield(opsNew,opsFields{iFN})
            opsNew.(opsFields{iFN}) = ops.(opsFields{iFN});
        else
            % If the substructure does exist, try looping through the
            % sub-structure, to check whether there are any missing bits
            try
                opsSubFields = fieldnames(ops.(opsFields{iFN}));
                
                for iSFN = 1:length(opsSubFields)
                    
                    if ~isfield(opsNew.(opssFields{iFN}),opsSubFields{iSFN})
                        
                        opsNew.(opsFields{iFN}).(opsSubFields{iSFN}) = ...
                            ops.(opsFields{iFN}).(opsSubFields{iSFN});
                    end
                end
            end
        end
        
    end
    
else
    opsNew = ops;
end
end

