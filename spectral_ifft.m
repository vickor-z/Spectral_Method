function [ f ] = spectral_ifft( hatf )
% һά
% �����׿ռ��ϵ���hatf�����ٸ���Ҷ��任���õ�����ռ��ֵf
N = (length(hatf) - 1) / 2;
f = fft(hatf);
j = 1 : 2 * N + 1;
f = f.*exp(2*N*(j-1)/(2*N+1)*pi*1i);

end

