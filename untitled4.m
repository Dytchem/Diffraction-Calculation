clear, clc, close all % test

l = 1 / 2;
L = 10;
N = 128;
Uc = @(x, y)1;
Ud = discretize(Uc, -l, l, -l, l, N, N);
Ud1 = myFFT2(Ud, -l, l, -l, l, -L, L, -L, L);
x = linspace(-L, L, N);
[x, y] = meshgrid(x, x);

surf(x, y, abs(Ud1))