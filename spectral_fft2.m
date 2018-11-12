function [hatu] = spectral_fft2(u)
%二维傅里叶变换
%获取物理空间u相应的谱hatu，用到了matlab自带的ifft2快速算法
%其中输入u为N*N的方阵，输出hatu也为N*N方阵，为u对应坐标的频率

%获取矩阵维度
N = size(u);
N = (N(1) - 1) / 2;

%对每个u的值我们需要先进行一次平移（x，y方向同时平移）
j = 0 : 2 * N;
u = u .* (exp(-2 * pi * N / (2 * N + 1) * j' * 1i) * exp(-2 * pi * N / (2 * N + 1) * j * 1i));

%快速变换
hatu = ifft2(u);
end

