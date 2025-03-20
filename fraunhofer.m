function [Uc1, Ud1, Ud] = fraunhofer(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N)
%FRAUNHOFER 计算连续分布在特定取样条件下的夫琅禾费衍射积分
%   小写字母为输入光场采样参数，大写字母为输出光场采样参数
%   lambda 波长
%   z 传播距离

k = 2 * pi / lambda;
Ud = discretize(Uc, xmin, xmax, ymin, ymax, m, n);
Ud1 = zeros(N, M);
x = linspace(xmin, xmax, m);
y = linspace(ymin, ymax, n);
[x, y] = meshgrid(x, y);
for i = 1:M
    X = (Xmax - Xmin) / (M - 1) * (i - 1) + Xmin;
    for j = 1:N
        Y = (Ymax - Ymin) / (N - 1) * (j - 1) + Ymin;
        Ud1(j, i) = sum(Ud.*exp(-1j*k/z*(X * x + Y * y)), "all");
    end
end

X = linspace(Xmin, Xmax, M);
Y = linspace(Ymin, Ymax, N);
[X, Y] = meshgrid(X, Y);
Ud1 = Ud1 * exp(1j*k*z) / (1j * lambda * z) .* exp(1j * k * (X.^2 + Y.^2) / (2 * z)) * (xmax - xmin) / (m - 1) * (ymax - ymin) / (n - 1);
Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);

end
