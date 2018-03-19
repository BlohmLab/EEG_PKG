function [ eeg, channels, filtT ] = PreProcessEEG( S, alignType, Fs, mu_rxn_samp )
%PRE Summary of this function goes here
%   Detailed explanation goes here
%Argument housekeeping
if nargin < 4
    mu_rxn_samp = 0;  
end


%Get fields from struct
s_fields = fields(S);
c_fields = fields(S.(s_fields{1}));

filtT = []; 
%Perform preprocessing operation
for stim = 1:length(s_fields)-1
    for cond = 1:length(c_fields)
        
        validBlock = unique(S.(s_fields{stim}).(c_fields{cond}).EEG(:,end));
        
        for b = 1:length(validBlock)
            
            %Grab current block
            currBlock = S.(s_fields{stim}).(c_fields{cond});
            data = currBlock.EEG(currBlock.EEG(:,end) == validBlock(b),:);
            data(:,2:6) = data(:,2:6)./sqrt(1000);
            Teeg = 1000*(data(end,1) - data(1,1))/length(data(:,1));
            
            %Find trials
            ind = find(data(:,7) == 3);
            trig = group(ind,1, 10);
            events = currBlock.T(:,:,validBlock(b))./1000;
            
            if ~xor(length(trig) ~= 100, length(trig) ~= 200)
                disp(['Problematic grouping!!! \n b = ' num2str(b) ' \n cond = ' num2str(cond) ' \n stim = ' num2str(stim)]);
                              
                %Get number of expected trials in this paradigm
                %100 add modifier if exceeds arm reach length
                maxTrials = 50 + (50 * (length(trig) > 100));
                
                %Get number of missing trials 
                numMissing = maxTrials - length(trig)/2;
                events(1:numMissing,:) = []; 
                
                t_append = [stim*ones(numMissing,1), cond*ones(numMissing,1),...
                    b*ones(numMissing,1) (2:numMissing+1)']; 
                filtT = [filtT; t_append]; 
                t_append = [];
                
            end
            
            %Get trial events, short fix EEG beginning after trial start,
            %delete it... 
            if (currBlock.T(1,1:size(currBlock.T,2),validBlock(b)) == 0 | ...
                    isnan(currBlock.T(1,1:size(currBlock.T,2),validBlock(b))))
                
                trig(1:2) = [];
                events(1,:) = [];
                
            end
            
           %% Filtering
            [nb, db] = butter(3, 2*Teeg./1000*[1 50]);
            [nn, dn] = butter(3, 2*Teeg./1000*[55 65], 'stop');
            for k = 1:5
                datafilt(:,k) = data(:,k+1);
                datafilt(:,k) = filtfilt(nn, dn, data(:,k+1));
                datafilt(:,k) = filtfilt(nb, db, datafilt(:,k));
                datafilt(:,k) = (datafilt(:,k) - mean(datafilt(:,k)))./std(datafilt(:,k));
                zdataf = abs(zscore(diff(data(:,k+1))));
                probch = find(sum(zdataf > 7.5) > 0);
            end
%             for i = 1:length(probch)
%                 disp(['Channel ' num2str(probch(i)-1) ' is problematic']);
%             end
            
            disp(['stim: ' num2str(stim) ' cond: ' num2str(cond) ' b: ' num2str(validBlock(b))])
            
           %% Event alignment

            %Align to selected event (-1 trial, since we ignore last one) 
            [block_eeg{b,cond,stim}, block_channels{b,cond,stim}] = ...
                AlignEEGData_v2([data(:,1) datafilt data(:,7)], trig,...
                events(:,2+alignType) - events(:,1), Fs, alignType, mu_rxn_samp);  
            
            datafilt = []; 
        end
    end
end

%Concatenate cell of EEG
% for i = 1:size(block_eeg,2) * size(block_eeg,3)
%     sub_eeg{1+ (i > size(block_eeg,2)), 2 - mod(i,2)} = vertcat([block_eeg{:,2 - mod(i,2), 1+ (i > size(block_eeg,2))}])'; 
%     sub_channels{1+ (i > size(block_eeg,2)), 2 - mod(i,2)} = vertcat([block_channels{:,2 - mod(i,2), 1+ (i > size(block_eeg,2))}])'; 
% end


%CHANGE FOR REMOVING LEARNING 
for stim = 1 : size(block_eeg,3)
    for cond = 1 : size(block_eeg,2) 
        sub_eeg{stim,cond} = vertcat([block_eeg{1:end,cond,stim}])'; 
        sub_channels{stim,cond} = vertcat([block_channels{1:end,cond,stim}])'; 
    end
end

eeg = sub_eeg; 
channels = sub_channels; 
end

