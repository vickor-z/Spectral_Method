function [u] = spectral_ifft2(hatu)
%二维逆傅里叶变换
%获取谱空间hatu相应物理空间的值u，用到了matlab自带的fft2快速算法
%其中输入hatu为N*N的方阵，输出u也为N*N方阵，为hatu对应坐标的值

%获取矩阵维度
N = size(hatu);
N = (N(1) - 1) / 2;

%先进行快速变换
u = fft2(hatu);

%变换后需要对每个值进行一次平移（x，y方向同时平移)
j = 0 : 2 * N;
u = u .* (exp(2 * pi * N / (2 * N + 1) * j' * 1i) * exp(2 * pi * N / (2 * N + 1) * j * 1i));
end

