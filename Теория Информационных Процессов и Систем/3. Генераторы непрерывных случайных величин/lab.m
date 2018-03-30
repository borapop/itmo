
a = 22695477;
c = 1;
m = 2^32;

X_size = 1000;

Y = zeros(1, X_size * 10);

for i = 2:size(Y, 2)
    Y(i) = mod(a * Y(i - 1) + c, m);
end

for i = 1:size(Y, 2)
    Y(i) = Y(i) / m;
end

syms x

F1 = @(x) (x >= 0) .* (x < 2) .* (0.125 .* x);
F2 = @(x) (x >= 2) .* (x < 4) .* (0.125 .* x - 0.25);
F3 = @(x) (x >= 4) .* (x < 6) .* (0.125 .* x - 0.5);
F4 = @(x) (x >= 6) .* (x <= 8) .* (0.125 .* x - 0.75);

F = @(x) F1(x) + F2(x) + F3(x) + F4(x);
F_max = 0.25;
a = 0;
b = 8;

X = zeros(1, X_size);
j = 1;
for i = 1:X_size
    eta = a + Y(j) * (b - a);
    while(Y(j + 1) * F_max > F(eta))
        j = j + 2;
        eta = a + Y(j) * (b - a);
    end
    X(i) = eta;
    j = j + 2;
end

%hist(X, 20);



M = integral(@(x) F(x) .* x, 0, 8)
D = integral(@(x) F(x) .* (x - M) .^ 2, 0, 8)

M_e = 0;
for i = 1:X_size
    M_e = M_e + X(i);
end
M_e = M_e / X_size

D_e = 0;
for i = 1:X_size
    D_e = D_e + (X(i) - M_e) ^ 2;
end
D_e = D_e / (X_size - 1)

alpha = 1 - 0.95;
u = 1.96; %U_0.975

radius = (D_e * u) / sqrt(X_size);
CI_a = M_e - radius
CI_b = M_e + radius

%correlation
C = zeros(1, 20);

for i = 1:20
    cov = 0;
    for j = 1:20 - i
        cov = cov + (X(j) - M_e) * (X(j + i) - M_e);
    end
    cov = cov / (20 - i);
    C(i) = cov / D_e;
end

x = linspace(0, 20, 20);
figure
plot(x, C, '-');

figure
plot(X(1:20 - 1), X(2:20), 'x');

figure
plot(X, linspace(0, 1, 1), '.');
fplot(F, [0, 8])
hold on
hist(X, 20)
hold off



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



