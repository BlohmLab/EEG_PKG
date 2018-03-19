function [ alignedEEG, channels ] = AlignEEGData_v2( data, trig, eventTime, Fs, alignType, mu_rxn_samp)
%ALIGNEEGDATA Summary of this function goes here
%   Inputs required:
% data - EEG data
%trial group - using trigger message system, indices of start and end of
%trial
% eventTime - time from trigger start to event required
%alignType - stimulus onset, or activity onset (0 and 1 respectively)
%mu_rxn_samp - median reaction time in samples for RT based events

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
    counter = 1; 
    %Ditch last trial due to high variability in end time
    for i = 1:length(trig)/2 - 1
        
        if ~isnan(alignEvent(i))
            if (i ~= length(trig)/2 -1)
                %This goes to just prior to the next trial (includes ITI) 
                curTrialInd = ind(trig(2*i-1)):ind(trig(2*i+1))-1; %:trig(2*i+1)-1);
            else
                curTrialInd = ind(trig(2*i-1)):size(data,1);
            end
            
            arr = curTrialInd(alignEvent(i)) - (250 + alignBaseOffset):curTrialInd(alignEvent(i)) + 500;
            eeg{k}(:,counter) = data(arr,k+1);
            
            counter = counter + 1;
        end
    end
end

alignedEEG = eeg{1}; 
channels = 0.25*(eeg{2} + eeg{3} + eeg{4} + eeg{5}); 

