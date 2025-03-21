function [Uc1, Ud1, Ud] = kirchhoff(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, Xmin, Xmax, Ymin, Ymax, M, N)
%KIRCHHOFF 计算连续分布在特定取样条件下的基尔霍夫衍射积分初步近似
%   小写字母为输入光场采样参数，大写字母为输出光场采样参数
%   lambda 波长
%   z 传播距离

k = 2 * pi / lambda;
[Ud, x, y] = discretize(Uc, xmin, xmax, ymin, ymax, m, n);
Ud1 = zeros(N, M);

for i = 1:M
    X = (Xmax - Xmin) / (M - 1) * (i - 1) + Xmin;
    for j = 1:N
        Y = (Ymax - Ymin) / (N - 1) * (j - 1) + Ymin;
        Ud1(j, i) = sum(Ud.*exp(1j*k*sqrt(z^2+(x - X).^2+(y - Y).^2)), "all");
    end
end

Ud1 = Ud1 / (1j * lambda * z) * (xmax - xmin) / (m - 1) * (ymax - ymin) / (n - 1);
Uc1 = interpolate(Ud1, Xmin, Xmax, Ymin, Ymax);

end
