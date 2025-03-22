%% 高斯拉盖尔光束（使角谱正确）
clear, clc, close all

lambda = 633e-9; % 波长
r = 5e-5; % 波片半径
w = 1e-5; % 参数
l = 5; % 参数
z = 0.0001; % 传播距离
Xmin = -2 * r;
Xmax = 2 * r;
Ymin = -2 * r;
Ymax = 2 * r; % 观察屏范围
nn = 1024; % 高分辨率

xmin = -r;
xmax = r;
ymin = -r;
ymax = r;
Uc = @(x, y)(x^2 + y^2 < r^2) .* ((x.^2 + y.^2) .* exp(-(x.^2 + y.^2)/w^2+1j*l*atan2(x, y)));


m = nn;
n = nn;
M = nn;
N = nn;

[~, Ud1, Ud] = fresnel_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
[~, Ud2] = fraunhofer_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
[~, Ud3] = jiaopu(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);

x = linspace(xmin, xmax, m);
y = linspace(ymin, ymax, n);
[x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);

figs = [];
work(x, y, Ud, "初始光场U_0(x_0,y_0)");
work(X, Y, Ud1, "菲涅尔衍射光场");
work(X, Y, Ud2, "夫琅禾费衍射光场");
work(X, Y, Ud3, "角谱传播计算光场");


function work(x, y, U, name)
figure
subplot(1, 2, 1)
surf(x, y, abs(U), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
xlabel("x")
ylabel("y")
zlabel("|U|")
title(name+" 振幅")
cmin = min(abs(U(:))); % 数据最小值
cmax = max(abs(U(:))); % 数据最大值
colorbar('Limits', [cmin, cmax]); % 固定颜色范围

subplot(1, 2, 2)
surf(x, y, mod(angle(U), 2*pi), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
xlabel("x")
ylabel("y")
zlabel("\phi")
title(name+" 相位")
cmin = 0; % 数据最小值
cmax = 2 * pi; % 数据最大值
colorbar('Limits', [cmin, cmax]); % 固定颜色范围
end