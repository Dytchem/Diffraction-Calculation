clear,clc,close all

m=2;n=10;r0=3e-10;E=6.4087e-19;

a=r0^((m-n));
alpha=E/((1-m/n)*a^(-m/(m-n)))
beta=m*alpha/(n*a)