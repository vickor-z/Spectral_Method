function [ f ] = spectral_ifft( hatf )
% 一维
% 对在谱空间上的谱hatf作快速傅里叶逆变换，得到物理空间的值f
N = (length(hatf) - 1) / 2;
f = fft(hatf);
j = 1 : 2 * N + 1;
f = f.*exp(2*N*(j-1)/(2*N+1)*pi*1i);

end

