function TestPreProcessEEG( S, alignType, Fs, mu_rxn_samp )
%PRE Summary of this function goes here
%   Detailed explanation goes here
%Argument housekeeping
if nargin < 4
    mu_rxn_samp = 0;
end

persistent sub;
global E T;

%Keep track of subject
if isempty(sub)
    sub = 1;
else
    sub = sub + 1;
end


s_file = {'an','ca'};
c_file = {'pr','pt'};
a_file = {'F'};


%Get fields from struct
s_fields = fields(S);
c_fields = fields(S.(s_fields{1}));

filtT = [];
%Perform preprocessing operation
for stim = 1:length(s_fields)-1
    for cond = 1:length(c_fields)
        
        validBlock = unique(S.(s_fields{stim}).(c_fields{cond}).EEG(:,end));
        
        for b = 1:length(validBlock)
            
            fileStr = [E(sub).subInitials{:} s_file{stim} c_file{cond} num2str(b) a_file{1} ...
                '_Dmat.mat'];
            
            %Grab current block
            currBlock = S.(s_fields{stim}).(c_fields{cond});
            data = currBlock.EEG(currBlock.EEG(:,end) == validBlock(b),:);
            data(:,2:6) = data(:,2:6)./sqrt(1000);
            Teeg = 1000*(data(end,1) - data(1,1))/length(data(:,1));
            
            %Find trials
            ind = find(data(:,7) == 3);
            trig = group(ind,1, 10);
            events = currBlock.T(:,:,b)./1000;
            
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
            
            disp(['stim: ' num2str(stim) ' cond: ' num2str(cond) ' b: ' num2str(b)])
            
            %% Event alignment
            
            %TESTING
            ConfirmAlignment([data(:,1) datafilt data(:,7)], trig,...
                events(:,2+alignType) - events(:,1), Fs, alignType, mu_rxn_samp,fileStr);
            
            datafilt = [];
        end
    end
end

end

