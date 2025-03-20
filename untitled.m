clear, clc, close all % test


U0 = @(x, y)x^2;


U = discretize(U0, -10, 10, -10, 10, 20, 10);


U1 = interpolate(U, -10, 10, -10, 10);

[X, Y] = meshgrid(linspace(-10, 10, 100), linspace(-10, 10, 100));
% X = X';
% Y = Y';


% surf(X, Y, U)
xlabel("x")
ylabel("y")
zlabel("U")
figure
surf(X, Y, U1(X, Y));