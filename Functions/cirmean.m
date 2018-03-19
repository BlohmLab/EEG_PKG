function [m, r, ci, sd] = cirmean(x, alp)
%
% x: sample angular values (vector) in rad (not deg)
% alp: confidence level for confidence interval (CI)
%
% m: mean angle
% r: mean resultant length
% ci: confidence limits on m
% sd: sample circular SD of m

if nargin<2, alp = 0.05; end

N = length(x);

S = sum(sin(x));
C = sum(cos(x));

r = sqrt(S.^2+C.^2)./N;
m = atan2(S, C);

if r < 0.53, 
    k = 2.*r + r.^3 + 5./6.*r.^5;
elseif (r >= 0.53 & r < 0.85),
    k = -.4 + 1.39.*r + .43.*(1-r);
elseif r >= 0.85,
    k = 1./(r.^3 - 4.*r.^2 + 3.*r);
end

lev = norminv((alp+1)/2,0,1);

ci(1) = m - lev./sqrt(N.*r.*k);
ci(2) = m + lev./sqrt(N.*r.*k);

sd = sqrt(2.*(1-r));