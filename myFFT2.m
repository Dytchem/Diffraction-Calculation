function Ud1 = myFFT2(Ud, xmin, xmax, ymin, ymax, Xmin, Xmax, Ymin, Ymax, M, N)
%MYFFT2 求离散分布Ud(X0,Y0)对应的连续分布Uc(x0,y0)的傅里叶变换Uc1(x,y)=∫∫Uc(x0,y0)exp(-i2pi(x0*x+y0*y))dx0dy0对应的离散分布Ud1(X,Y)
%   xmin,xmax,ymin,ymax 原函数的积分范围
%   Xmin,Xmax,Ymin,Ymax 像函数的计算范围
%   M,N 像函数的网格数量

m = size(Ud, 2);
n = size(Ud, 1);
kx = (xmax - xmin) / (m - 1);
ky = (ymax - ymin) / (n - 1);
kX = (Xmax - Xmin) / (M - 1);
kY = (Ymax - Ymin) / (N - 1);
ix = 0:m - 1;
iy = 0:n - 1;
Ix = 0:M - 1;
Iy = 0:N - 1;
IX = [Ix, -m:-1];
IY = [Iy, -n:-1];
[ix, iy] = meshgrid(ix, iy);
[Ix, Iy] = meshgrid(Ix, Iy);
[IX, IY] = meshgrid(IX, IY);

U1 = Ud .* exp(-1j*pi*(kx * kX * ix.^2 + ky * kY * iy.^2 + 2 * kx * Xmin * ix + 2 * ky * Ymin * iy));
U1(n+N, m+M) = 0;
U2 = exp(1j*pi*(kx * kX * IX.^2 + ky * kY * IY.^2));
Ud1 = ifft2(fft2(U1).*fft2(U2)); % 卷积定理
% Ud1 = cconv2(U1, U2);

Ud1 = Ud1(1:N, 1:M) .* exp(-1j*pi*(kx * kX * Ix.^2 + ky * kY * Iy.^2 + 2 * kX * xmin * Ix + 2 * kY * ymin * Iy + 2 * xmin * Xmin + 2 * ymin * Ymin)) * kx * ky;

end
