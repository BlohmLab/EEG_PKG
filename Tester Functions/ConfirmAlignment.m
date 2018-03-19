function ConfirmAlignment( data, trig, eventTime, Fs, alignType, mu_rxn_samp,fileStr)
%ALIGNEEGDATA Summary of this function goes here
%   Inputs required:
% data - EEG data
%trial group - using trigger message system, indices of start and end of
%trial
% eventTime - time from trigger start to event required
%alignType - stimulus onset, or activity onset (0 and 1 respectively)
%mu_rxn_samp - median reaction time in samples for RT based events

persistent timeError;
persistent badStrings;
persistent zeroTrials;
global outString totalError totalZeros;


if isempty(badStrings)
    badStrings = {};
end

if isempty(zeroTrials)
    zeroTrials = 0;
end

%Load in the D matrix first element
load(fileStr,'D');

numTrials = length(trig)/2;
trialsSkipped = 100 - numTrials;
validD = D(trialsSkipped + 1 : end);

%Find expected trial lengths for all included trials
trialDiff = cellfun(@(x) x.tITI - x.tFixation,validD,'un',false);
trialLength = cell2mat(trialDiff)./1000;


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

badAlign = false;
throwError = false;


%Ditch last trial due to high variability in end time
for i = 1:length(trig)/2 - 1
    
    if (~isnan(alignEvent(i)))
        
        %This goes to just prior to the next trial (includes ITI)
        triggerInd = ind(trig(2*i-1)):ind(trig(2*i));
        
        
        %TESTING FUNCTIONALITY
        % Test # 1: Check trial length match
        expectedLength = trialLength(i);
        currentLength = data(triggerInd(end),1) - data(triggerInd(1),1);
        delta = expectedLength - currentLength;
        
        %Convert delta into samples
        f_delta = delta/512;
        
        if (f_delta > 3)
            throwError = true;
        end
        
        %Routine B: Check error in sample and build up
        timeError = [timeError; delta];
        
        
        %Routine C: For when SPon Erraneously appears as 0.
        if (eventTime(i) == 0)
            throwError = true;
            
            %Routine D: To confirm the number of total zeros for match
            %checking
            zeroTrials = zeroTrials + 1;
            disp(zeroTrials);
            
        end
        
        
        
    end
end

if (badAlign)
    badStrings = [badStrings; fileStr];
end

if throwError
    badStrings = [badStrings; fileStr ' Zero'];
end

outString = badStrings;
totalError = timeError;
totalZeros = zeroTrials;
end





