clear, clc, close all % test

Uc = @(x, y)x.^0;
lambda = 1;
z = 1;
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
m = 32;
n = 32;

[Uc1, Ud1] = kirchhoff(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, xmin, xmax, ymin, ymax, m, n);
[Uc2, Ud2] = fresnel(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, xmin, xmax, ymin, ymax, m, n);
[Uc3, Ud3] = fraunhofer(Uc, xmin, xmax, ymin, ymax, m, n, lambda, z, xmin, xmax, ymin, ymax, m, n);

X = linspace(xmin, xmax, m);
Y = linspace(ymin, ymax, n);
[X, Y] = meshgrid(X, Y);

figs = [];
figs(end+1) = subplot(2, 2, 1);
surf(X, Y, abs(Uc(X, Y)));
title("原光场")
figs(end+1) = subplot(2, 2, 2);
surf(X, Y, abs(Ud1));
title("基尔霍夫衍射")
figs(end+1) = subplot(2, 2, 3);
surf(X, Y, abs(Ud2));
title("菲涅尔衍射")
figs(end+1) = subplot(2, 2, 4);
surf(X, Y, abs(Ud3));
title("夫琅禾费衍射")
for f = figs
    xlabel(f, "x")
    ylabel(f, "y")
    zlabel(f, "U")
end
