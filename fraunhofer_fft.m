function [Uc1, Ud1, Ud] = fraunhofer_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N)
%FRAUNHOFER_FFT 计算连续分布在特定取样条件下的夫琅禾费衍射积分（FFT版）
%   小写字母为输入光场采样参数，大写字母为输出光场采样参数
%   lambda 波长
%   z 传播距离

k = 2 * pi / lambda;
Ud = discretize(Uc, xmin, xmax, ymin, ymax, m, n);
% x = linspace(xmin, xmax, m);
% y = linspace(ymin, ymax, n);
% [x, y] = meshgrid(x, y);
X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);

t = 1 / (lambda * z); % 缩放因子
Ud1 = myFFT2(Ud, xmin, xmax, ymin, ymax, Xmin*t, Xmax*t, Ymin*t, Ymax*t, M, N);

Ud1 = Ud1 * exp(1j*k*z) / (1j * lambda * z) .* exp(1j*k/(2 * z)*(X.^2 + Y.^2));
Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);

end
