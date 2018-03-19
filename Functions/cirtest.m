function p = cirtest(x1, x2, method)
%
% performs bi-variate statistical tests on data with
% hypothesized von Mises distribution (equivalent to
% standard distribution for non-circular data)
%
% x1: sample 1
% x2: sample 2
% method: 'Watson-Williams' --> tests if mean(x1) different mean(x2)
%
% p: p-level of difference between samples


if method(1) == 'W',
    
    N1 = length(x1);
    N2 = length(x2);
    N = N1 + N2;
    
    S1 = sum(sin(x1));
    C1 = sum(cos(x1));
    S2 = sum(sin(x2));
    C2 = sum(cos(x2));
    S = S1 + S2;
    C = C1 + C2;
    
    R1 = sqrt(S1.^2+C1.^2);
    R2 = sqrt(S2.^2+C2.^2);
    r1 = R1/N1;
    r2 = R2/N2;
    R = sqrt(C.^2 + S.^2);
    r = (R1 + R2)/N;

    if r < 0.53,
        k = 2.*r + r.^3 + 5./6.*r.^5;
    elseif (r >= 0.53 & r < 0.85),
        k = -.4 + 1.39.*r + .43.*(1-r);
    elseif r >= 0.85,
        k = 1./(r.^3 - 4.*r.^2 + 3.*r);
    end

    if k > 2,
        F = (1 + 3/8/k)*((R1 + R2 - R)/(N - R1 - R2));
        p = 1 - fcdf(F, 1, N-2);
    else
        disp('The concentration parameter of the von Mises distribution is < 2. Please chose another method.')
    end

end