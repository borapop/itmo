syms x

Fa = 0.3 .* x;
Fb = 0.2 + 0 .* x;

M = int(x .* Fa, 0, 2) + int(x .* Fb, 2, 4)
D = int(Fa .* (x - M) .^ 2, 0, 2) + int(Fb .* (x - M) .^ 2, 0, 4)

%{
%M[X]
M = 0;
for i = 1:7
    M = M + X(i) * P(i);
end
disp 'M[X]';
disp(M);

%D[X]
D = 0;
for i = 1:7
    D = D + P(i) * (X(i) - M) ^ 2;
end
disp 'D[X]';
disp(D);

E = 2.45 * sqrt(D / 7);
disp 'CI'
disp(M - E)
disp(M + E)


%correlation
corr = 1:7;

for i = 1:7
    cov = 0;
    for j = 1:7-i
        cov = cov + (X(j) - M) * (X(j + i) - M);
    end
    cov = cov / 7;
    corr(i) = cov / D;
end

figure
plot(X, corr, '-');

figure;
plot(X, P, '-');

figure;
x = linspace(0, 6, 6);
plot(x, X(1:6), 'o', x, X(2:7), 'o')
%}



