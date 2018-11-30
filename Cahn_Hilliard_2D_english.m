%Cahn-Hilliard equation
%parameter initialization
%2N + 1 is the number of the points, [a,b] is the domain, dt is the time step, eps is the parameter, u0=u(x,0)
%t is the beginning time, T is the ending time
NN = [20 ,40, 80, 160];
tt = [0.25,0.1,0.05,0.02];
for Np = 1:4
N = 100;
a = 0; b = 2 * pi ; dt = tt(Np); h = (b - a) / (2 * N + 1);
x = h * (0 : 2 * N)'; y = h * (0 : 2 * N)';
eps = 0.01; u0 = 0.05 * sin(x) * sin(y');
%u0 = 2 * rand(2 * N + 1) - 1;
t = 0; T = 2;

%initialization
%w is the nonlinear term, hatu is the spectral of u, A is the (p^2 + q^2) matrix
%W1,W2 is the translation matrix for fft
wp = u0 - u0.^3;
j = -N : N;
j = j .^ 2;
A = ones(2 * N + 1,1) * j + j' * ones(1, 2 * N + 1);
up = u0;
j = 0 : 2 * N;
W1 = (exp(-2 * pi * N / (2 * N + 1) * j' * 1i) * exp(-2 * pi * N / (2 * N + 1) * j * 1i));
W2 = (exp(2 * pi * N / (2 * N + 1) * j' * 1i) * exp(2 * pi * N / (2 * N + 1) * j * 1i));
hatup = ifft2(u0 .* W1);
flag = 0;

%Iteration, *p shows the value of time t, *n shows the value of time t+dt
while t < T 
    
    
    %calculate the spectral of w
    hatwp = ifft2(wp .* W1);
    
    %update the spectral of u
    hatun = (hatup  + dt * A .* hatwp ) ./ (1 + eps * dt * A.^2);
    
    %transform the spectral to physical space, update u and w
    un = real(fft2(hatup) .* W2);
    wn = un - un.^3;
    
    %record residual
    res = max(abs(un - up));
    
    %time developing
    t = t + dt; up = un; wp = wn; hatup = hatun;
    
    %plot images
    %if floor(t * 7.5) == flag
    %    subplot(4,4,flag + 1);
    %    imagesc(x,y,up);
    %    title(['t = ',num2str(t)]);
    %    flag = flag + 1;
    %end
end
   subplot(2,2,Np);
   contour(x,y,up);
   title(['dt = ',num2str(dt), ', N = ', num2str(N)]);
end
%take in the end point
X=[x;b];
Y=[y;b];
U=[up up(1:2*N+1,1);up(1,1:2*N+1) up(1,1)];

%plot contour line
%contour(X,Y,U);


        