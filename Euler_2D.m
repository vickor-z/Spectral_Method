%test1 谱方法在基上对任意N是精确的
N = 10;
a = 0;
b = 2*pi;
x = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
y = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
df = 2*sin(x')*sin(y);


%取f的谱
hatf = spectral_fft2(df);

%计算u的谱
j = -N:N;
j = j.^2;
tr = ones(2*N + 1, 2*N + 1);
for k = 1:2*N +1
    tr(k,1:2*N+1) = j(1:2*N+1)+ j(k);
end
tr(N+1,N+1)=1;
hatu = hatf./tr;
hatu(N+1,N+1) = 0;

%注意到，ifft是先平移后变换，fft是先变换后平移
u=spectral_ifft2(hatu);


%作图
uu = sin(x')*sin(y);
[X,Y]=meshgrid(x,y);
res = max(abs(u-uu));
mesh(X,Y,uu);

        