%解Allen-Cahn方程
%参数初始化
%2N + 1是点的个数，[a,b]是计算区间，dt是时间步长，x是该周期内的前2N+1个点，y为相应的y轴坐标，Lx为周期长度，eps是参数，u0=u(x,0)
%t为初始时刻，T为结束时刻
N = 20;
a = 0; b = 2 * pi ; dt = 1 / (2 * N + 1) / (2 * N + 1);
x = a : (b - a) / (2 * N + 1) : b - (b - a) / (2 * N + 1);
y = a : (b - a) / (2 * N + 1) : b - (b - a) / (2 * N + 1);
Lx = b - a; eps = 0.01; u0 = sin(x') * sin(y);

t = 0; T = 1;

%初始化，其中w为非线性项的值，hatu为u的谱，A为迭代过程中的方阵，但仅需点乘即可
w = u0 - u0.^3;
hatu = spectral_fft2(u0);
j = -N : N;
j = j .^ 2;
A = zeros(2 * N + 1, 2 * N + 1);
for k = 1 : 2 * N + 1
    A(1 : 2 * N + 1, k) = 1 ./ (1 + dt * eps * (j(1 : 2 * N + 1) + j(k)));
end
u = u0;

%迭代
while t < T
    %为了记录残量变化，我们需要记录上一个时刻的u值
    up = u;
    
    %计算非线性项w的谱
    hatw = spectral_fft2(w);
    
    %更新u的谱
    hatu = (hatu  + dt.*hatw ).*A;
    
    %将hatu投影回物理空间，更新u的值和w的值，在物理空间中，我们仅取实部即可
    u = real(spectral_ifft2(hatu));
    w = u - u.^3;
    
    %时间发展
    t = t + dt;
    
    %记录残量，这里定义residual为u^(n+1)-u^n的L无穷范数（而不是u-u^n的L无穷范数，所以并不是数学意义上的残量）
    res = max(abs(u - up));
    
    %绘制动态等高线图
    %contour(x,y,u);
    %pause(0.01);
end

%绘制等高线图
contour(x,y,u);
title({'The solution of Allen Cahn Equation with fourier galerkin method';'when T=1,epsilon=0.01 and u0=sinxsiny'});


        