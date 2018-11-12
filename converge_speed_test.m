%test2 测试收敛速度，验证res=e^-N
res=zeros(100,1);
for N = 1:100;
a = 0;
b = 2*pi;
x = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
df = f(x);
j = 1:2*N+1
df=df.*exp(-2*N*(j-1)/(2*N+1)*pi*i);
hatf = ifft(df);


n = (j-N-1).*(j-N-1);
n(N+1) = 1;
hatu = hatf./n;
hatu(N+1) = 1;


u=fft(hatu);
u=u.*exp(2*N*(j-1)/(2*N+1)*pi*i);

y=3./(5-4.*cos(x));
res(N)=max(abs(y-u));
end
%plot(x,y,x,u);
%注意要取log10(res)
semilogy(1:100,(res));
        