%% 快速傅里叶变换积分（FFT）、角谱传播公式的对比
clear, clc, close all

r = 1e-5;
Uc = @(x, y)x.^2 + y.^2 <= r^2;
lambda = 633e-9;
z = 0.001;
xmin = -r;
xmax = r;
ymin = -r;
ymax = r;
Xmin = -r * 10;
Xmax = r * 10;
Ymin = -r * 10;
Ymax = r * 10;
nn = 1024; % 超高分辨率
m = nn;
n = nn;
M = nn;
N = nn;

ts = [];
tic
[Uc1, Ud1, Ud] = fresnel_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc; tic
[Uc2, Ud2] = fraunhofer_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc; tic
[Uc3, Ud3] = jiaopu(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
ts(end+1) = toc;

x = linspace(xmin, xmax, m);
y = linspace(ymin, ymax, n);
[x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);

figure
figs = [];
figs(end+1) = subplot(2, 2, 1);
surf(x, y, abs(Ud), "EdgeColor", "none");
title("初始光场U_0(x_0,y_0)")
figs(end+1) = subplot(2, 2, 2);
surf(X, Y, abs(Ud1), "EdgeColor", "none");
title("菲涅尔衍射（FFT）："+ts(1)+"秒")
figs(end+1) = subplot(2, 2, 3);
surf(X, Y, abs(Ud2), "EdgeColor", "none");
title("夫琅禾费衍射（FFT）："+ts(2)+"秒")
figs(end+1) = subplot(2, 2, 4);
surf(X, Y, abs(Ud3), "EdgeColor", "none");
title("角谱传播公式（FFT）："+ts(3)+"秒")

for f = figs
    xlabel(f, "x")
    ylabel(f, "y")
    zlabel(f, "|U|")
end
