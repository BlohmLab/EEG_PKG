function [ alignedEEG, channels ] = AlignEEGData( data, trig, eventTime, Fs, alignType, mu_rxn_samp)
%ALIGNEEGDATA Summary of this function goes here
%   Inputs required:
% data - EEG data
%trial group - using trigger message system, indices of start and end of
%trial
% eventTime - time from trigger start to event required
%alignType - stimulus onset, or activity onset (0 and 1 respectively)
%mu_rxn_samp - median reaction time in samples for RT based events
%opt - 'laplacian' for surface laplacian

alignTargEvent = eventTime;

ind = find(data(:,end) == 3); 
switch(alignType)
    case 0
        alignEvent = round(alignTargEvent*Fs);
        alignBaseOffset = 0;
    case 1
        alignEvent = round(alignTargEvent*Fs);
        alignBaseOffset = round(mu_rxn_samp);
end

for k = 1:5
    for i = 1:length(trig)/2 - 1
        curTrialInd = ind(trig(2*i-1):trig(2*i));
        arr = curTrialInd(alignEvent(i)) - (250 + alignBaseOffset):curTrialInd(alignEvent(i)) + 500;
        eeg{k}(:,i) = data(arr,k+1);
    end
end

alignedEEG = eeg{1}; 
channels = 0.25*(eeg{2} + eeg{3} + eeg{4} + eeg{5}); 

