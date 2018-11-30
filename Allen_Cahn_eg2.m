%解Allen-Cahn方程
%参数初始化
%2N + 1是点的个数，[a,b]是计算区间，dt是时间步长，x是该周期内的前2N+1个点，y为相应的y轴坐标，Lx为周期长度，eps是参数，u0=u(x,0)
%t为初始时刻，T为结束时刻
N = 20;
a = 0; b = 2 * pi ; dt = 0.001; h = (b - a) / (2 * N + 1);
x = h * (0 : 2 * N)'; y = h * (0 : 2 * N)';
eps = 0.01; %u0 = sin(x) * sin(y');
u0 = 2 * rand(2 * N + 1) - 1;
t = 0; T = 1;

%初始化，其中w为非线性项的值，hatup为u在t时刻的谱，A为迭代过程中的方阵，但仅需点乘即可
%W1为快速变换的平移矩阵, W2为快速逆变换的平移矩阵
wp = u0 - u0.^3;
hatup = spectral_fft2(u0);
j = -N : N;
j = j .^ 2;
A = zeros(2 * N + 1, 2 * N + 1);
for k = 1 : 2 * N + 1
    A(1 : 2 * N + 1, k) = 1 ./ (1 + dt * eps * (j(1 : 2 * N + 1) + j(k)));
end
up = u0;
j = 0 : 2 * N;
W1 = (exp(-2 * pi * N / (2 * N + 1) * j' * 1i) * exp(-2 * pi * N / (2 * N + 1) * j * 1i));
W2 = (exp(2 * pi * N / (2 * N + 1) * j' * 1i) * exp(2 * pi * N / (2 * N + 1) * j * 1i));

%迭代，迭代过程中，*p表示t时刻的值，*n表示t+dt时刻的值
while t < T 
    %计算非线性项w的谱
    hatwp = ifft2(wp .* W1);
    
    %更新u的谱
    hatun = (hatup  + dt .* hatwp ).*A;
    
    %将hatu投影回物理空间，更新u的值和w的值，在物理空间中，我们仅取实部即可
    un = real(fft2(hatup) .* W2);
    wn = un - un.^3;
    
    %记录残量，这里定义residual为u^(n+1)-u^n的L无穷范数（而不是u-u^n的L无穷范数，所以并不是数学意义上的残量）
    res = max(abs(un - up));
    
    %时间发展
    t = t + dt; up = un; wp = wn; hatup = hatun;
    
    X=[x;b];
    Y=[y;b];
    U=[up up(1:2*N+1,1);up(1,1:2*N+1) up(1,1)];

    %绘制动态等高线图
    contour(X,Y,U);
    %pause(0.001);
end

%扩展坐标，将周期另一端加入
X=[x;b];
Y=[y;b];
U=[up up(1:2*N+1,1);up(1,1:2*N+1) up(1,1)];

%绘制等高线图
contour(X,Y,U);


        