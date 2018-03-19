function [R, p] = cirregress(x, y, method)
%
% calculates correlation coefficient between two circular variables
%
% x, y: samples (rad)
% method: 'Uniform' --> for uniformely distributed angles
%         'Ranks' --> rank correlation coefficient (for N > 8)
%         'Linear' --> x = angular variable, y = linear variable
%         'Parametric' --> alternative to the Uniform (using chi-square
%         test for large N)
%
% R: correlation coefficient

if length(y) ~= length(x),
    disp('X and Y must have the same length')
else
    if method(1) == 'U',
        N = length(x);

        dp = y - x;
        dn = y + x;

        rp = sum(cos(dp - mean(dp)))/N;
        rn = sum(cos(dn - mean(dn)))/N;
        
        R = max(rp, rn);
        p = 1 - raylcdf(N*R^2, 1);
    elseif method(1) == 'R',
        N = length(x);
        e = 2*pi/N;
        
        [s, i] = sort(y);
        
        xt = (1:N)*e;
        yt = i*e;
        
        dp = yt - xt;
        dn = yt + xt;
        
        rp = sum(cos(dp - mean(dp)))/N;
        rn = sum(cos(dn - mean(dn)))/N;

        R = max(rp, rn);
        p = 1 - (1/(1 - exp(-(N-1)*R^2)))^2;
    elseif method(1) == 'P',
        N = length(x);
        rcc = circorr(cos(x), cos(y));
        rcs = circorr(cos(x), sin(y));
        rsc = circorr(sin(x), cos(y));
        rss = circorr(sin(x), sin(y));
        r1 = circorr(cos(x), sin(x));
        r2 = circorr(cos(y), sin(y));
        
        R = sqrt((rcc^2+rcs^2+rsc^2+rss^2 + 2*(rcc*rss+rcs*rsc)*r1*r2 -...
            2*(rcc*rcs+rsc*rss)*r2 - 2*(rcc*rsc+rcs*rss)*r1)/((1-r1^2)*(1-r2^2)));
        p = 1 - chi2cdf(N*R^2, 4);
    elseif method(1) == 'L',
        N = length(x);
        ryc = circorr(y, cos(x));
        rys = circorr(y, sin(x));
        rcs = circorr(cos(x), sin(x));
        R = sqrt((ryc^2+rys^2-2*ryc*rys*rcs)/(1-rcs^2));
        p = 1 - chi2cdf(N*R^2, 2);
    end
end
    
    