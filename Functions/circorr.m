function r = circorr(x1, x2),

N = length(x1);
c = sum((x1 - mean(x1)).*(x2 - mean(x2)))/(N-1);
s1 = sqrt(sum(x1 - mean(x1))^2)/(N-1);
s2 = sqrt(sum(x2 - mean(x2))^2)/(N-1);

r = c/(s1*s2);