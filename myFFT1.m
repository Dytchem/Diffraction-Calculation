function Ud1 = myFFT1(Ud, xmin, xmax, Xmin, Xmax, M)
%MYFFT1 求离散分布Ud(X0)对应的连续分布Uc(x0)的傅里叶变换Uc1(x)=∫Uc(x0)exp(-i2pi*x0*x)dx0对应的离散分布Ud1(X)

m = length(Ud);
kx = (xmax - xmin) / (m - 1);
kX = (Xmax - Xmin) / (M - 1);
ix = 0:m - 1;
Ix = 0:M - 1;
IX = [Ix, -m:-1];

U1 = Ud .* exp(-1j*pi*(kx * kX * ix.^2 + 2 * kx * Xmin * ix));
U1(1, m+M) = 0;
U2 = exp(1j*pi*kx*kX.*IX.^2);

Ud1 = ifft(fft(U1).*fft(U2)); % 卷积定理
% Ud1 = cconv(U1, U2, m+M);

Ud1 = Ud1(1:M) .* exp(-1j*pi*(kx * kX * Ix.^2 + 2 * kX * xmin * Ix + 2 * xmin * Xmin)) * kx;

end
