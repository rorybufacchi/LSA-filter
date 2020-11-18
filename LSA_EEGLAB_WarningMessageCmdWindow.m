function [] = LSA_EEGLAB_WarningMessageCmdWindow()
%[] = LSA_EEGLAB_WarningMessageCmdWindow()
%   Displays a warning message about the interpretability of the results in
%   the command window

warningText = {'IMPORTANT: interpretation of results',...
     ' ',...
     'Not all of the assumptions that make LSA interpretable will hold at each timepoint', ...
     'Guaranteed interpretability only exists for the indicated channels. Therefore, this function will output ', ...
     'three EEG structures, and save three separate EEG datasets:', ...
     '  - 1) EEG_LSA: This dataset contains the full filtered timecourse. ', ...
     '        IMPORTANTLY, there is no guarantee that each timepoint is interpretable. ', ...
     '  - 2) EEG_LSA_nan: This dataset contains only the interpretable parts of the filtered timecourse. ', ...
     '        All timepoints for which the assumptions do not hold are set to NaN. ', ...
     '        Note that you might not be able to plot this output in EEGLAB due to its handling of NaNs. ', ...
     '  - 3) EEG_LSA_selection: This dataset contains the average of the data within the time window selected by the user. ', ...
     '        Hence, this dataset contains only 1 timepoint. All channels for which the assumptions do not hold are set to NaN. Therefore, the ', ...
     '        output data are fully interpretable. Note that plotting this output in EEGLAB becomes difficult as there is only 1 timepoint.' }

for iL = 1:length(warningText)
    disp(warningText{iL})
end


end

