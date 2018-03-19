function [TFR, PLOC] = makeTFR2(S,freqVec,timeVec,Fs,width,TimeInt, ver)
% D. Cheyne, September 2006
%
% function TFR = makeTFR(trialData,freqVec,timeVec,Fs,width,TimeInt);
%
%
% based on old routine 4D toolbox routine from Ole Jensen to get TFR but 
% includes baseline adjustment in one call.  
% function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width);
%
% Calculates the average of a time-frequency energy representation of
% multiple trials using a Morlet wavelet method.                            
%
% Input
% -----
% S    : signals = time x Trials  
% timeVec    : absolute time values for trial        
% freqVec    : frequencies over which to calculate TF energy        
% Fs   : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
% TimeInt: array specifying start and end times for baseline, e.g. [-4 -3]
%
% Output
% ------
% t    : time
% f    : frequency
% B    : phase-locking factor = frequency x time
%     
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------
%

% Edited by Andrew Hamilton, to produce text-based progress messages.

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%======================================================================
% Modified by Gunnar Blohm to compute z-scores of data re. baseline
%======================================================================

if nargin < 7
    ver = 1; 
end
S = S';
LNR50 = 0; 

%timeVec = (1:size(S,2))/Fs;  
tidx = find(timeVec >= TimeInt(1) & timeVec <= TimeInt(2));
B = zeros(length(freqVec),size(S,2)); 
pow = zeros(length(freqVec),size(S,2)); 
disp('TFR Status: Beginning calculation ...')

for ii=1:size(S,1)          
    %fprintf(1,'%d ',ii); 
    for jj=1:length(freqVec)
        
        % get power and phase for this bin
        [t,z] = energyvec(freqVec(jj),detrend(S(ii,:)),Fs,width);
        b = mean(t(tidx));
        s = std(t(tidx));
        B(jj,:) = (t-b) + B(jj,:);
        PLOC(jj,:,ii) = z; 
        SS(jj,ii) = s;
    end
	% text output of calculation progress
	if mod(ii,10) == 0,
		disp(['TFR Status: ' sprintf('%0.0f', ii/size(S,1)*100 ) '%']);
	end
end
switch ver
    case 1
        TFR = B/size(S,1)./repmat(mean(SS,2),1,length(B(1,:)));
    case 2
         TFR = B/size(S,1);%./repmat(mean(SS,2),1,length(B(1,:)));
end
disp(['TFR Status: Finished']);


function [y,z] = energyvec(f,s,Fs,width)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
%
% 

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
y = conv(s,m);
z0 = y./abs(y); %Normalize complex morlet (not necessary?) - doesnt make a difference to angle
z1 = angle(z0); %Compute phase angle
 
%Subtract out the base frequency, or expected frequency cycling... 



y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2)); %%Remove padding and get signal
z = z1(ceil(length(m)/2):length(z1)-floor(length(m)/2)); 



function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = 1/(st*sqrt(2*pi));

% yes, this has been corrected to the right scaling across frequencies...
% Complex Morlet Wavelet, can project real and complex domains (real and imag functions)
y = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t); 
