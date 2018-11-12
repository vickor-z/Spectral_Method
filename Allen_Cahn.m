%解Allen-Cahn方程
%参数初始化
%2N + 1是点的个数，[a,b]是计算区间，dt是时间步长，x是该周期内的前2N+1个点，Lx为周期长度，eps是参数，u0=u(x,0)
%t为初始时刻，T为结束时刻
N = 20;
a = 0; b = 2 * pi; dt = 1 / (2 * N + 1) / (2 * N + 1);
x = a : (b - a) / (2 * N + 1) : b - (b - a) / (2 * N + 1);
Lx = b - a; eps = 0.01; u0 = sin(x - pi);

t = 0; T = 1;

%初始化，其中w为非线性项的值，hatu为u的谱，A为迭代过程中的对角阵，我们在此作向量处理
w = u0 - u0.^3;
hatu = spectral_fft(u0);
j = -N : N;
A = 1 ./ (1 + eps .* j .* j .* dt);
u = u0;

%迭代
while t<T
    %为了记录残量变化，我们需要记录上一个时刻的u值
    up = u;
    
    %计算非线性项w的谱
    hatw = spectral_fft(w);
    
    %更新u的谱
    hatu = (hatu  + dt.*hatw ).*A;
    
    %将hatu投影回物理空间，更新u的值和w的值，在物理空间中，我们仅取实部即可
    u = real(spectral_ifft(hatu));
    w = u - u.^3;
    
    %时间发展
    t=t+dt;
    
    %记录残量，这里定义residual为u^(n+1)-u^n的L无穷范数（而不是u-u^n的L无穷范数，所以并不是数学意义上的残量）
    res = max(abs(u - up));
    
end

plot(x,u);




        