s = tf('s');

a12 = 0.1;
a21 = 0.05;
V = 1;

A = [-a12-0.01, a21;
    a12, -a21-0.01];
B = [1; 0];
C = [1, 0];

c = [0; 0.0];
x= 0.0;

nt = 50;
N = 10000;
t = linspace(0, nt, N);
dt = nt/(N);

r = ones(1, N)*0.1;
y = zeros(1, N+1);
c1 = zeros(1, N);

k =22;


for i=1:N
    u = k*((r(i)-y(i))>0)*(r(i)-y(i));
    c = c + dt * ( A*c + B*u );

    x = x + dt * 0.1*( -x + c(1) );
    % ym = Clin*x;
    c1(i) = C*c;
    y(i+1) = interp1(t, c1, t(i)-0.1, 'linear', 0);
end

figure
plot(t, y(2:end))


%%

G = ss(A, B, C, 0);
H = exp(-s*0.1)*tf(1, [10, 1]);

figure
nyquist(G*H);
% hold on
% nyquist(k/2);
