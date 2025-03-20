clear, clc, close all % test

a = rand(1, 10);
b = rand(1, 10);

c = cconv(a, b, 10);
d = ifft(fft(a).*fft(b));

disp([c; d])