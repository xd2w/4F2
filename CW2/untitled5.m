close all;

s = tf('s');

m = 0.001;
cp = 2;
cv = 0.1;

A1 = [
    -cv/m,    -(cp+0.1)/m;
      1, 0;
    ];

A = [
    -cv/m,    -cp/m;
      1, 0;
    ];

A2 = [
    -cv/m,    -(cp-0.1)/m;
      1, 0;
    ];


B = [1/m; 0];
C = [1, 30];
D= 0;
[n, ~] = size(A);

lambda = 20;
G_lin = ss(A, B, C, D);

Alin = A;
Blin = B;
Clin = C;

ep = 1e-6
cvx_begin sdp
    cvx_solver SDPT3

    variable P(n, n) symmetric
    variable a
    variable b
    variable ep
    % minimise(ep)
    % minimise(b-a)
    % LMI1 = P > 0;
    LMI2 = a < 0;
    LMI3 = ep > 0;
    LMI4 = [A1'*P +  P*A1 + 2*lambda*P - a*(C'*C) + ep*eye(n), P*B-C';
            B'*P-C, -b] < 0;
    LMI5 = [A2'*P +  P*A2 + 2*lambda*P- a*(C'*C) + ep*eye(n), P*B-C';
            B'*P-C, -b] < 0;
cvx_end

display(['excess pasivity of' num2str(a)])

% none of the poles P must be negative
eig(P)


% 0 poles left of lambda=20 should work anyway
eig(A1)
eig(A2)

a_lin = a
b_lin = b

% nt = 10;
% N = 1000000;
% t = linspace(0, nt, N);
% dt = nt/(N);
% 
% q = ones(1, N)*0;
% rho =10;
% 
% ym = zeros(1, N); 
% y = zeros(1, N+1);
% 
% z = [0; 0; 0.00];
% x = [0; 0];
% 
% for i=1:N
%     x = x + dt * ( Alin*x + Blin*(y(i)+q(i) ) );
%     ym(i) = Clin*x;
% end

%%

A = [
    -10,    0,  10;
      0, -100, 100;
      0,    0,  -1;
    ];

B = [0; 0; 1];
C = [0, 0, -1];
D= 0;
[n, ~] = size(A);

G_nonlin = ss(A, B, C, D);

Anonlin = A;
Bnonlin = B;

lambda = 20;
ep = 1e-6;

cvx_begin sdp
    cvx_solver SDPT3

    variable Y(n, n) symmetric
    variable Z(1, n)
    variable a
    variable b
    % variable ep
    % maximise(-b)
    % minimise(b-a)
    % LMI1 = Y>0
    LMI2 = b >= 0;
    % LMI4 = Y*A' +  A*Y + Z'*B' + B*Z + 2*lambda*Y + ep*eye(n) < 0
    % LMI5 = Y*A' +  A*Y + 2*lambda*Y + ep*eye(n) < 0
    LMI6 = Y*A' +  A*Y + Z'*B' + B*Z + ep*eye(n) < 0;
    LMI4 = [[Y*A' +  A*Y + Z'*B' + B*Z + 2*lambda*Y + ep*eye(n), B-Y*C', Y*C'];
                [ B'-C*Y, -b , D'];
                [ C*Y , D , a]] < 0;
    LMI5 = [[Y*A' +  A*Y + 2*lambda*Y + ep*eye(n), B-Y*C', Y*C'];
                [ B'-C*Y, -b , D'];
                [ C*Y , D , a]] < 0;
cvx_end

ep

eig(Y)

b_nonlin = b
a_nonlin = 1/a

K = Z/Y

rho = 1;

nt = 10;
N = 10000;
t = linspace(0, nt, N);
dt = nt/(N);

q = ones(1, N)*0;

r = ones(1, N)*0; 
y = zeros(1, N+1);

z = [0; 0; 0.01];
x = [0; 0];

for i=1:N
    z = z + dt * ( G_nonlin.A*z + G_nonlin.B*(tanh(K*z) + r(i)) );
    y(i) = -z(3);
end

figure
plot(t, y(2:end))
% hold on
% plot(t, r(1:end))

%%

nt = 10;
N = 1000000;
t = linspace(0, nt, N);
dt = nt/(N);

q = ones(1, N)*.8;
rho =0.1;

ym = zeros(1, N); 
y = zeros(1, N+1);

z = [0; 0; 0.001];
x = [0; 0];

for i=1:N
    x = x + dt * ( Alin*x + Blin*(y(i)+q(i) ) );
    ym(i) = Clin*x;
    z = z + dt * ( Anonlin*z + Bnonlin*(tanh(K*z) + rho*ym(i)) );
    y(i+1) = -z(3);
end

figure
% plot(t, y(2:end))
% hold on
plot(t, ym(1:end))
% ylim([0, 1])
% xlim([0, 0.1])

a_lin
b_nonlin

b_lin 
a_nonlin