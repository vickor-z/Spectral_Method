%test1 �׷����ڻ��϶�����N�Ǿ�ȷ��
N = 10;
a = 0;
b = 2*pi;
x = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
y = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
df = 2*sin(x')*sin(y);


%ȡf����
hatf = spectral_fft2(df);

%����u����
j = -N:N;
j = j.^2;
tr = ones(2*N + 1, 2*N + 1);
for k = 1:2*N +1
    tr(k,1:2*N+1) = j(1:2*N+1)+ j(k);
end
tr(N+1,N+1)=1;
hatu = hatf./tr;
hatu(N+1,N+1) = 0;

%ע�⵽��ifft����ƽ�ƺ�任��fft���ȱ任��ƽ��
u=spectral_ifft2(hatu);


%��ͼ
uu = sin(x')*sin(y);
[X,Y]=meshgrid(x,y);
res = max(abs(u-uu));
mesh(X,Y,uu);

        