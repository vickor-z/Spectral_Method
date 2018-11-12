function [u] = spectral_ifft2(hatu)
%��ά�渵��Ҷ�任
%��ȡ�׿ռ�hatu��Ӧ����ռ��ֵu���õ���matlab�Դ���fft2�����㷨
%��������hatuΪN*N�ķ������uҲΪN*N����Ϊhatu��Ӧ�����ֵ

%��ȡ����ά��
N = size(hatu);
N = (N(1) - 1) / 2;

%�Ƚ��п��ٱ任
u = fft2(hatu);

%�任����Ҫ��ÿ��ֵ����һ��ƽ�ƣ�x��y����ͬʱƽ��)
j = 0 : 2 * N;
u = u .* (exp(2 * pi * N / (2 * N + 1) * j' * 1i) * exp(2 * pi * N / (2 * N + 1) * j * 1i));
end

