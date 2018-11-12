function [ hatf ] = spectral_fft( f )
% 对在离散点上的函数f做快速傅里叶变换，得到hatf
N = (length(f) - 1) / 2;
j = 1 : 2 * N + 1;
f = f .* exp( -2 * N * (j - 1) / (2 * N + 1) * pi * 1i);
hatf = ifft(f);
end

