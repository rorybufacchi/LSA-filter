function [] = eegplugin_LSAfilt(fig, try_strings, catch_strings)


% find main menu
toolsmenu = findobj(fig, 'tag', 'tools');

% build command for menu callback
cmd = [ '[EEG_LSA com ops EEG_LSA_nan EEG_LSA_selection] = pop_LSAfilt_EEGLAB(EEG); ' ...
        '[ALLEEG EEG CURRENTSET]                         = eeg_store(ALLEEG, EEG_LSA); ' ...
        '[ALLEEG EEG CURRENTSET]                         = eeg_store(ALLEEG, EEG_LSA_nan); ' ...
        '[ALLEEG EEG CURRENTSET]                         = eeg_store(ALLEEG, EEG_LSA_selection); '];
 
finalcmd = [ try_strings.no_check cmd ]; 
finalcmd = [ finalcmd 'LASTCOM = ''' cmd ''';' ]; 
finalcmd = [ finalcmd catch_strings.store_and_hist ]; 
 
% add new submenu
submenu = uimenu( toolsmenu, 'label', 'LSA filter', 'callback', finalcmd);