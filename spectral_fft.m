function [ hatf ] = spectral_fft( f )
% ������ɢ���ϵĺ���f�����ٸ���Ҷ�任���õ�hatf
N = (length(f) - 1) / 2;
j = 1 : 2 * N + 1;
f = f .* exp( -2 * N * (j - 1) / (2 * N + 1) * pi * 1i);
hatf = ifft(f);
end

