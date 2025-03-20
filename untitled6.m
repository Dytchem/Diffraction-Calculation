clear, clc, close all % test

l = 1 / 2;
L = 10;
N = 10000;
Uc = @(x)1 + 0 .* x;
Ud = Uc(linspace(-l, l, N));
Ud1 = myFFT1(Ud, -l, l, -L, L, N);
x = linspace(-L, L, N);

plot(x, abs(Ud1))