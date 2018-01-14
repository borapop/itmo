a = 22695477;
c = 1;
m = 2^32;
X0 = 1;
X(1) = 0;
X(1000) = 0;

Y(1) = 100;
Y(1000) = 0;

for i = 1:1000
    Y(i + 1) = mod(a * Y(i) + c, m);
    X(i) = mod(a * Y(i) + c, m) / 2^32;
    disp(X(i))
end

%M[X]

M = 0;
P = 1 / 1000;
for i = 1:1000
    M = M + X(i) * P;
end
disp 'M[X]';
disp(M);

%D[X]
D = 0;
for i = 1:1000
    D = D + P * (X(i) - M) ^ 2;
end
disp 'D[X]';
disp(D);

%k, ro

k = 1:999;
po = 1:999;



corr = 1:20;

for i = 1:20
    cov = 0;
    for j = 1:20-i
        cov = cov + (X(j) - M) * (X(j + i) - M);
    end
    cov = cov / (20 - i);
    corr(i) = cov / D;
end
x = linspace(0, 20, 20);
plot(x, corr, '-');
figure;
plot(X(1:999), X(2:1000), '.');



figure;
h = histogram(X);
h.NumBins = 10;
h.Normalization = 'probability';