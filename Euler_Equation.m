%test1 �׷����ڻ��϶�����N�Ǿ�ȷ��
N = 10;
a = 0;
b = 2*pi;
x = a:(b-a)/(2*N+1):b-(b-a)/(2*N+1);
df = 4*sin(2*x);
j=1:2*N+1
%ƽ��f
df=df.*exp(-2*N*(j-1)/(2*N+1)*pi*i);
%ȡf����
hatf = ifft(df);

%����u����
n = (j-N-1).*(j-N-1);
n(N+1) = 1;
hatu = hatf./n;
hatu(N+1) = 0;

%ע�⵽��ifft����ƽ�ƺ�任��fft���ȱ任��ƽ��
u=fft(hatu);
u=u.*exp(2*N*(j-1)/(2*N+1)*pi*i);

%��ͼ
y=sin(2*x);
plot(x,y,x,u);
res=max(abs(y-u))



        