close all

a1 = 22695477;
c1 = 1;
m1 = 2^32;

a2 = 1103515245;
m2 = 2^31;
c2 = 12345;

Y_size = 1000;

Y1 = zeros(1, Y_size);
Y1(1) = 1;
Y2 = zeros(1, Y_size);
Y2(1) = 1;

for i = 2:Y_size
    Y1(i) = mod(a1 * Y1(i - 1) + c1, m1);
end

for i = 1:Y_size
    Y1(i) = Y1(i) / m1;
end

for i = 2:Y_size
    Y2(i) = mod(a2 * Y2(i - 1) + c2, m2);
end

for i = 1:Y_size
    Y2(i) = Y2(i) / m2;
end

X_size = Y_size;

X1 = zeros(1, X_size);
X2 = zeros(1, X_size);

M = 1;
D = 0.5;

for i = 1:X_size
    X1(i) = sqrt(-2 * log( Y1(i) ) ) * cos(2 * pi * Y2(i));
    X1(i) = M + sqrt(D) * X1(i);
    
    X2(i) = sqrt(-2 * log( Y1(i) ) ) * sin(2 * pi * Y2(i));
    X2(i) = M + sqrt(D) * X2(i);
end


M1_e = count_M_e(X1, X_size)
D1_e = count_D_e(M1_e, X1, X_size)

M2_e = count_M_e(X2, X_size)
D2_e = count_D_e(M2_e, X2, X_size)

CI1 = count_CI(M1_e, D1_e, X_size)
CI2 = count_CI(M2_e, D2_e, X_size)

plot_correlation(M1_e, D1_e, X1, 1000)
plot_correlation(M2_e, D2_e, X2, 1000)

figure
histogram(X1)
figure
histogram(X2)

function plot_correlation(M_e, D_e, X, size)
    C = zeros(1, size);

    for i = 1:size
        cov = 0;
        for j = 1:size - i
            cov = cov + (X(j) - M_e) * (X(j + i) - M_e);
        end
        cov = cov / (size - i);
        C(i) = cov / D_e;
    end
    x = linspace(0, size, size);
    figure
    plot(x, C, '-');

    figure
    plot(X(1:size - 1), X(2:size), '.');
end

function CI = count_CI(M_e, D_e, X_size)
    u = 1.96; %U_0.975
    radius = (D_e * u) / sqrt(X_size);
    CI_a = M_e - radius;
    CI_b = M_e + radius;
    CI = [CI_a, CI_b];
end


function M_e = count_M_e(X, X_size)
    M_e = 0;
    for i = 1:X_size
        M_e = M_e + X(i);
    end
    M_e = M_e / X_size;
end

function D_e = count_D_e(M_e, X, X_size)
    D_e = 0;
    for i = 1:X_size
        D_e = D_e + (X(i) - M_e) ^ 2;
    end
    D_e = D_e / (X_size - 1);
end

%{

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
%}



