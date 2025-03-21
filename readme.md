# 相干光场在自由空间中传播衍射的数值计算

[![View Diffraction-Calculation on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/180450-diffraction-calculation)

---

## 三种衍射积分形式

- 基尔霍夫衍射初步近似（简单离散积分）

$$
U\left( x,y \right) =\frac{1}{\mathrm{j}\lambda z}\iint_{-\infty}^{+\infty}{U_0\left( x_0,y_0 \right) \mathrm{e}^{\mathrm{j}k\sqrt{z^2+\left( x-x_0 \right) ^2+\left( y-y_0 \right) ^2}}\mathrm{d}x_0\mathrm{d}y_0}
$$

- 菲涅尔衍射（简单离散积分 or 快速傅里叶变换快速卷积）

$$
U\left( x,y \right) =\frac{\mathrm{e}^{\mathrm{j}kz}}{\mathrm{j}\lambda z}\iint_{-\infty}^{+\infty}{U_0\left( x_0,y_0 \right) \mathrm{e}^{\mathrm{j}k\frac{\left( x-x_0 \right) ^2+\left( y-y_0 \right) ^2}{2z}}\mathrm{d}x_0\mathrm{d}y_0}
$$

- 夫琅禾费衍射（简单离散积分 or 快速傅里叶变换快速卷积）

$$
U\left( x,y \right) =\frac{\mathrm{e}^{\mathrm{j}kz}}{\mathrm{j}\lambda z}e^{\mathrm{j}\frac{k}{2z}\left( x^2+y^2 \right)}\iint_{-\infty}^{+\infty}{U_0\left( x_0,y_0 \right) \mathrm{e}^{-\mathrm{j}\frac{k}{z}\left( x_0x+y_0y \right)}\mathrm{d}x_0\mathrm{d}y_0}
$$

---

## 两种积分运算方式

- 简单离散积分

时间复杂度$O(N^2)$

- 快速傅里叶变换快速卷积

![积分离散化](src/f1.png)

其中

![快速傅里叶变换和快速卷积](src/f2.png)

时间复杂度$Nlog(N)$

代码请见 [myFFT1.m](myFFT1.m) 和 [myFFT2.m](myFFT2.m)

---

## 效率对比

- [32*32](task1.m)

![32*32](src/task1.png)

- [128*128](task2_1.m)

![128*128](src/task2_1.png)

- [1024*1024](task2_2.m)
![1024*1024](src/task2_2.png)

---

## 衍射样例

- [双缝干涉](task3_1.m)

![双缝干涉初始光场](src/task3_1.png)
![双缝干涉衍射光场](src/task3_2.png)
![双缝干涉衍射光场](src/task3_3.png)]

- [振幅型余弦光栅](task4.m)

![振幅型余弦光栅初始光场](src/task4_1.png)
![振幅型余弦光栅衍射光场](src/task4_2.png)
![振幅型余弦光栅衍射光场](src/task4_3.png)

- [泊松亮斑](task5.m)

![泊松亮斑初始光场](src/task5_1.png)
![泊松亮斑衍射光场](src/task5_2.png)
![泊松亮斑衍射光场](src/task5_3.png)

- [高斯拉盖尔光束](task6.m)

![高斯拉盖尔光束初始光场](src/task6_1.png)
![高斯拉盖尔光束衍射光场](src/task6_2.png)
![高斯拉盖尔光束衍射光场](src/task6_3.png)
![高斯拉盖尔光束衍射光场](src/task6_4.png)
