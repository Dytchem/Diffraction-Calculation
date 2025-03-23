%% 验证圆孔的夫琅禾费衍射公式
clear, clc, close all

%% 模型参量
lambda = 500e-9;
z = 1;
r = 1e-3;
U0 = 1;
nn = 1024;
xmin = -r;
xmax = -xmin;
ymin = -r;
ymax = -ymin;
Xmin = xmin;
Xmax = -Xmin;
Ymin = ymin;
Ymax = -Ymin;
x = linspace(xmin, xmax, nn);
y = linspace(ymin, ymax, nn);
[x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, nn);
Y = linspace(Ymin, Ymax, nn);
[X, Y] = meshgrid(X, Y);
k = 2 * pi / lambda;

%% 理论解
U1 = @(x, y)exp(1j*k*(z + (x.^2 + y.^2) / (2 * z))) * r ./ (1j * sqrt(x.^2+y.^2)) .* besselj(1, 2*pi*r*sqrt(x.^2+y.^2)/(lambda * z));
U1 = U1(X, Y);
work(X, Y, U1, "理论解");

%% 数值解
[~, U2] = fraunhofer_fft(@(x, y)x.^2+y.^2 <= r^2, xmin, xmax, ymin, ymax, nn, nn, lambda, z, Xmin, Xmax, Ymin, Ymax, nn, nn);
work(X, Y, U2, "数值解");

%% 误差函数
work(X, Y, U2-U1, "数值-理论");

%% 绘图函数
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