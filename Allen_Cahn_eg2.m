%��Allen-Cahn����
%������ʼ��
%2N + 1�ǵ�ĸ�����[a,b]�Ǽ������䣬dt��ʱ�䲽����x�Ǹ������ڵ�ǰ2N+1���㣬yΪ��Ӧ��y�����꣬LxΪ���ڳ��ȣ�eps�ǲ�����u0=u(x,0)
%tΪ��ʼʱ�̣�TΪ����ʱ��
N = 20;
a = 0; b = 2 * pi ; dt = 1 / (2 * N + 1) / (2 * N + 1);
x = a : (b - a) / (2 * N + 1) : b - (b - a) / (2 * N + 1);
y = a : (b - a) / (2 * N + 1) : b - (b - a) / (2 * N + 1);
Lx = b - a; eps = 0.01; u0 = sin(x') * sin(y);

t = 0; T = 1;

%��ʼ��������wΪ���������ֵ��hatuΪu���ף�AΪ���������еķ��󣬵������˼���
w = u0 - u0.^3;
hatu = spectral_fft2(u0);
j = -N : N;
j = j .^ 2;
A = zeros(2 * N + 1, 2 * N + 1);
for k = 1 : 2 * N + 1
    A(1 : 2 * N + 1, k) = 1 ./ (1 + dt * eps * (j(1 : 2 * N + 1) + j(k)));
end
u = u0;

%����
while t < T
    %Ϊ�˼�¼�����仯��������Ҫ��¼��һ��ʱ�̵�uֵ
    up = u;
    
    %�����������w����
    hatw = spectral_fft2(w);
    
    %����u����
    hatu = (hatu  + dt.*hatw ).*A;
    
    %��hatuͶӰ������ռ䣬����u��ֵ��w��ֵ��������ռ��У����ǽ�ȡʵ������
    u = real(spectral_ifft2(hatu));
    w = u - u.^3;
    
    %ʱ�䷢չ
    t = t + dt;
    
    %��¼���������ﶨ��residualΪu^(n+1)-u^n��L�������������u-u^n��L����������Բ�������ѧ�����ϵĲ�����
    res = max(abs(u - up));
    
    %���ƶ�̬�ȸ���ͼ
    %contour(x,y,u);
    %pause(0.01);
end

%���Ƶȸ���ͼ
contour(x,y,u);
title({'The solution of Allen Cahn Equation with fourier galerkin method';'when T=1,epsilon=0.01 and u0=sinxsiny'});


        