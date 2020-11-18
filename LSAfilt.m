function [flt, lbd] = LSAfilt(eeg, ref, mtd)
% LSAfilt Local Adaptive spatial filtering.
%
% USAGE
%   [FLT, LBD] = LSAfilt(EEG, REF, MTD)
%
% INPUT
%   EEG - Point-by-Trials-by-Channels EEG matrix
%   REF - Can be either a Point-by-Trial signal matrix or a channel index
%         of the EEG matrix (in that case EEG(:,:,REF) is used as reference
%         signal)
%   MTD - Method option: 'standard' | 'adaptive'

[nPnt, nTrl, nChn] = size(eeg);

if isscalar(ref)
    ref = eeg(:,:,ref);
end
ref = ref(:,:,ones(nChn,1));

refCpy = ref;
refAve = mean(refCpy, 2);
refCpy = refCpy - refAve;

eegCpy = eeg;
eegCpyAve = mean(eegCpy, 2);
eegCpy = eegCpy - eegCpyAve;

switch mtd

    case 'adaptive'
        lbd = sum(refCpy.*eegCpy, 2) ./ sum(refCpy.*refCpy, 2);

    case 'standard'
        lbd = ones(nPnt, 1, nChn);

    otherwise
        error('Case not found!');

end

flt = eeg - lbd(:,ones(nTrl,1),:) .* ref;
end