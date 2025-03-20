%% 直接求和积分（SUM）与快速傅里叶变换积分（FFT）的对比
clear, clc, close all

Uc = @(x, y)x^2 + y^2 <= (1e-3)^2;
lambda = 633e-9;
z = 100;
xmin = -1e-3;
xmax = 1e-3;
ymin = -1e-3;
ymax = 1e-3;
Xmin = -0.1;
Xmax = 0.1;
Ymin = -0.1;
Ymax = 0.1;
nn = 32; % 低分辨率
m = nn;
n = nn;
M = nn;
N = nn;

ts = [];
tic
[Uc1, Ud1, Ud] = kirchhoff(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc; tic
[Uc2, Ud2] = fresnel(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc; tic
[Uc3, Ud3] = fresnel_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc; tic
[Uc4, Ud4] = fraunhofer(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc; tic
[Uc5, Ud5] = fraunhofer_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc;

x = linspace(xmin, xmax, m);
y = linspace(ymin, ymax, n);
[x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);

figs = [];
figs(end+1) = subplot(3, 2, 1);
surf(x, y, abs(Ud));
title("初始光场U_0(x_0,y_0)")
figs(end+1) = subplot(3, 2, 2);
surf(X, Y, abs(Ud1));
title("基尔霍夫衍射（SUM），用时："+ts(1)+"秒")
figs(end+1) = subplot(3, 2, 3);
surf(X, Y, abs(Ud2));
title("菲涅尔衍射（SUM），用时："+ts(2)+"秒")
figs(end+1) = subplot(3, 2, 4);
surf(X, Y, abs(Ud3));
title("菲涅尔衍射（FFT），用时："+ts(3)+"秒")
figs(end+1) = subplot(3, 2, 5);
surf(X, Y, abs(Ud4));
title("夫琅禾费衍射（SUM），用时："+ts(4)+"秒")
figs(end+1) = subplot(3, 2, 6);
surf(X, Y, abs(Ud5));
title("夫琅禾费衍射（FFT），用时："+ts(5)+"秒")

for f = figs
    xlabel(f, "x")
    ylabel(f, "y")
    zlabel(f, "|U|")
end
