function [Uc1, Ud1, Ud] = fresnel_fft(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N)
%FRESNEL_FFT 计算连续分布在特定取样条件下的菲涅尔衍射积分（FFT版）
%   小写字母为输入光场采样参数，大写字母为输出光场采样参数
%   lambda 波长
%   z 传播距离

k = 2 * pi / lambda;
[Ud, x, y] = discretize(Uc, xmin, xmax, ymin, ymax, m, n);

X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);

t = 1 / (lambda * z); % 缩放因子
Ud1 = myFFT2(Ud.*exp(1j*k/(2 * z)*(x.^2 + y.^2)), xmin, xmax, ymin, ymax, Xmin*t, Xmax*t, Ymin*t, Ymax*t, M, N);

Ud1 = Ud1 .* exp(1j*k*z+1j*k/(2 * z)*(X.^2 + Y.^2)) / (1j * lambda * z);
Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);

end
