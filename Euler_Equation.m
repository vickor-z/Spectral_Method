%test1 谱方法在基上对任意N是精确的
N = 10;
a = 0;
b = 2*pi;
x = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
df = 4*sin(2*x);
j=1:2*N+1
%平移f
df=df.*exp(-2*N*(j-1)/(2*N+1)*pi*i);
%取f的谱
hatf = ifft(df);

%计算u的谱
n = (j-N-1).*(j-N-1);
n(N+1) = 1;
hatu = hatf./n;
hatu(N+1) = 0;

%注意到，ifft是先平移后变换，fft是先变换后平移
u=fft(hatu);
u=u.*exp(2*N*(j-1)/(2*N+1)*pi*i);

%作图
y=sin(2*x);
plot(x,y,x,u);
res=max(abs(y-u))



        