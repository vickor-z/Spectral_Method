%Allen-Cahn equation
%parameter initialization
%2N + 1 is the number of the points, [a,b] is the domain, dt is the time step, eps is the parameter, u0=u(x,0)
%t is the beginning time, T is the ending time
N = 20;
a = 0; b = 2 * pi ; dt = 0.01; h = (b - a) / (2 * N + 1);
x = h * (0 : 2 * N)'; y = h * (0 : 2 * N)';
eps = 0.01; %u0 = sin(x) * sin(y');
u0 = 2 * rand(2 * N + 1) - 1;
t = 0; T = 10;

%initialization
%w is the nonlinear term, hatu is the spectral of u, A is the iteration matrix
%W1,W2 is the translation matrix for fft
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

%Iteration, *p shows the value of time t, *n shows the value of time t+dt
while t < T 
    %calculate the spectral of w
    hatwp = ifft2(wp .* W1);
    
    %update the spectral of u
    hatun = (hatup  + dt .* hatwp ).*A;
    
    %transform the spectral to physical space, update u and w
    un = real(fft2(hatup) .* W2);
    wn = un - un.^3;
    
    %record residual
    res = max(abs(un - up));
    
    %time developing
    t = t + dt; up = un; wp = wn; hatup = hatun;
    
    %plot images
    contour(x,y,up);
    pause(0.001);
end

%take in the end point
X=[x;b];
Y=[y;b];
U=[up up(1:2*N+1,1);up(1,1:2*N+1) up(1,1)];

%plot contour line
contour(X,Y,U);


        