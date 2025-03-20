%% 计算双缝干涉条纹（叠加法计算）
clear, clc, close all

lambda = 633e-9; % 波长
a = 1e-5; % 缝宽
d = 1e-4; % 缝间距
c = 1; % 缝长度
z = 1; % 传播距离
Xmin = -0.01;
Xmax = 0.01;
Ymin = -1;
Ymax = 1; % 观察屏范围
% nn = 128; % 高分辨率

m = 512;
n = 512;
M = 16;
N = 2000;

% 左缝
xmin = -a / 2 - d / 2;
xmax = a / 2 - d / 2;
ymin = -c / 2;
ymax = c / 2;
Uc = @(x, y)ones(size(x));
[~, Ud1_1] = fresnel_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
[~, Ud2_1] = fraunhofer_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);

% 右缝
xmin = -a / 2 + d / 2;
xmax = a / 2 + d / 2;
ymin = -c / 2;
ymax = c / 2;
Uc = @(x, y)ones(size(x));
[~, Ud1_2] = fresnel_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);
[~, Ud2_2] = fraunhofer_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N);

Ud1 = Ud1_1 + Ud1_2;
Ud2 = Ud2_1 + Ud2_2;

% 原光场
xmin = -a / 2 - d / 2;
xmax = a / 2 + d / 2;
ymin = -c / 2;
ymax = c / 2;
Uc = @(x, y)abs(x) >= (d - a) / 2;
Ud = discretize(Uc, xmin, xmax, ymin, ymax, m, n);


x = linspace(xmin, xmax, m);
y = linspace(ymin, ymax, n);
[x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);

figs = [];
work(x, y, Ud, "初始光场");
work(X, Y, Ud1, "菲涅尔衍射光场");
work(X, Y, Ud2, "夫琅禾费衍射光场");


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
