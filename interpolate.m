function Uc = interpolate(Ud, xmin, xmax, ymin, ymax)
%INTERPOLATE 对给定离散分布插值成连续分布函数
%   Uc(x,y) 连续分布函数 Ud(X,Y) 离散分布矩阵
%   xmin,xmax,ymin,ymax x,y方向分布范围

M = size(Ud, 2);
N = size(Ud, 1);
X = linspace(xmin, xmax, M);
Y = linspace(ymin, ymax, N);
[X, Y] = meshgrid(X, Y);

Uc = @(x, y)interp2(X, Y, Ud, x, y, "spline");

end
