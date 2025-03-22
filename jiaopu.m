function [Uc1, Ud1, Ud] = jiaopu(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N, F, MM, NN)
%JIAOPU 计算连续分布在特定取样条件下的衍射积分（角谱传播）
%   小写字母为输入光场采样参数，大写字母为输出光场采样参数
%   lambda 波长
%   z 传播距离
%   角谱最大空间频率范围，默认取1/lambda
%   MM,NN 角谱的离散点数，默认取为M,N
if nargin <= 15
    F = 1 / lambda;
end
if nargin <= 16
    MM = M;
    NN = N;
end

k = 2 * pi / lambda;
Ud = discretize(Uc, xmin, xmax, ymin, ymax, m, n);

A = myFFT2(Ud, xmin, xmax, ymin, ymax, -F, F, -F, F, MM, NN); % 转化为频域量
XX = linspace(-F, F, MM);
YY = linspace(-F, F, NN);
[XX, YY] = meshgrid(XX, YY);
A = A .* exp(1j*k*z*sqrt(1-(lambda * XX).^2-(lambda * YY).^2)); % 角谱传播
surf(XX, YY, abs(A), "EdgeColor", "none")
xlabel("f_x")
ylabel("f_y")
zlabel("A")
title("角谱图")

Ud1 = myFFT2(A, -F, F, -F, F, -Xmax, -Xmin, -Ymax, -Ymin, M, N); % 转化回空域量
Ud1 = Ud1(end:-1:1, end:-1:1);
Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);

end
