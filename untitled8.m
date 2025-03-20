clear,clc,close all

x=[1,2];
y=[1,2,3];
[x,y]=meshgrid(x,y);

U=@(x,y)x*10+y;
disp(U(x,y))