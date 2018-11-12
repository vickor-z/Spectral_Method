function [hatu] = spectral_fft2(u)
%��ά����Ҷ�任
%��ȡ����ռ�u��Ӧ����hatu���õ���matlab�Դ���ifft2�����㷨
%��������uΪN*N�ķ������hatuҲΪN*N����Ϊu��Ӧ�����Ƶ��

%��ȡ����ά��
N = size(u);
N = (N(1) - 1) / 2;

%��ÿ��u��ֵ������Ҫ�Ƚ���һ��ƽ�ƣ�x��y����ͬʱƽ�ƣ�
j = 0 : 2 * N;
u = u .* (exp(-2 * pi * N / (2 * N + 1) * j' * 1i) * exp(-2 * pi * N / (2 * N + 1) * j * 1i));

%���ٱ任
hatu = ifft2(u);
end

