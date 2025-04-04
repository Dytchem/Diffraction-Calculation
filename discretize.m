function [Ud, X, Y] = discretize(Uc, xmin, xmax, ymin, ymax, M, N)
%DISCRETIZE 对给定连续分布离散取样
%   Uc(x,y) 连续分布函数 Ud(X,Y) 离散分布矩阵
%   xmin,xmax,ymin,ymax x,y方向分布范围
%   M,N x,y方向离散取点数

% 未矢量化的方法
% Ud = zeros(N, M);
% for i = 1:M
%     Ud(:, i) = arrayfun(@(y)Uc((xmax - xmin)/(M - 1)*(i - 1)+xmin, y), linspace(ymin, ymax, N));
% end

% 默认已矢量化
X = linspace(xmin, xmax, M);
Y = linspace(ymin, ymax, N);
[X, Y] = meshgrid(X, Y);

Ud = Uc(X, Y);

end