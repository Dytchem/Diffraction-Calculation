clear,clc,close all

A=[1,2,5];
B=[1,2,3];

cconv(A,B,3)
ifft(fft(A).*fft(B))