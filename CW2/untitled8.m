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

% Au = [
%     -10,    0,  10;
%       0, -100, 100;
%       0,    0,  -1;
%     ];

B = [1/m; 0];
C = [1, 30];
D= 0;
[n, ~] = size(A);

lambda = 20;
G_lin = ss(A, B, C, D);

Alin = A;
Blin = B;
Clin = C;

% nyquist(G)
% 
% figure
% rlocus(G)
ep = 1e-6
cvx_begin sdp
    cvx_solver SDPT3

    variable P(n, n) symmetric
    variable a
    variable b
    % variable ep
    minimise(b-a)
    % maximise(-ep)
    % LMI1 = P > 0;
    LMI2 = a < 0;
    LMI3 = ep > 0;
    % LMI4 = [A1'*P +  P*A1 - a*eye(n) + ep*eye(n), P*B-C';
    %         B'*P-C, -b] < 0;
    % LMI5 = [A2'*P +  P*A2 - a*eye(n) + ep*eye(n), P*B-C';
    %         B'*P-C, -b] < 0; 2*lambda*P
    % LMI6 = A'*P +  P*A + ep*eye(n)< 0;
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

G = ss(A, B, C, D);

% nt = 1;
% N = 10000;
% t = linspace(0, nt, N);
% dt = nt/(N);
% 
% q = ones(1, N)*0.1;
% 
% rho = 1;
% 
% 
% r = zeros(1, N); 
% y = zeros(1, N+1);
% 
% z = [0; 0; 0];
% x = [0; 0];
% 
% for i=1:N
%     x = x + dt * ( Alin*x + Blin*(q(i) ) );
%     ym = Clin*x;
%     r(i+1) = ym;
% end
% 
% figure
% plot(t, y(2:end))
% hold on
% plot(t, r(2:end))


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

ep = 1e-6;
lambda = 20;

A = [-10 0 10; 0 -100 100; 0 0 -1];
B = [0; 0; 1];
C = [0 0 -1];
D = 0;

cvx_begin sdp
    cvx_solver SDPT3

    variable Y(3,3) symmetric
    variable Z(1,3)
    variable a
    variable b
    % minimise(b-a)
    % minimise(ep)
    LMI1 = [Y*A' + A*Y + 2*lambda*Y + ep*eye(3) + B*Z + Z'*B', B - Y*C', Y*C';
            B' - C*Y, -b, 0;
            C*Y, 0, a] <= 0;
    LMI2 = [Y*A' + A*Y + 2*lambda*Y + ep*eye(3), B - Y*C', Y*C';
            B' - C*Y, -b, 0;
            C*Y, 0, a] <= 0;
    LMI3 = Y*A' + A*Y + Z'*B' + B*Z+ep*eye(3) <= 0;
    LMI4 = b > 0

cvx_end

K = Z/Y;

ep

eig(Y)

b_nonlin = b
a_nonlin = 1/a


rho = 1;

% nt = 10;
% N = 10000;
% t = linspace(0, nt, N);
% dt = nt/(N);
% 
% q = ones(1, N);
% 
% r = ones(1, N)*0.1; 
% y = zeros(1, N+1);
% 
% z = [0; 0; 0];
% x = [0; 0];
% 
% for i=1:N
%     z = z + dt * ( G_nonlin.A*z + G_nonlin.B*(tanh(K*z) + r(i)) );
%     y(i) = -z(3);
% end
% 
% figure
% plot(t, y(2:end))
% hold on
% plot(t, r(1:end))

%%

nt = 10;
N = 10000;
t = linspace(0, nt, N);
dt = nt/(N);

q = ones(1, N)*1;
q(N/2:end) = 0;
% q = 2*sin(t);

rho = 0.1;


r = zeros(1, N); 
y = zeros(1, N+1);

z = [0; 0; 0];
x = [0; 0.1];

for i=1:N
    z = z + dt * ( Anonlin*z + Bnonlin*(tanh(K*z) + r(i)) );
    y(i) = -z(3);
    x = x + dt * ( Alin*x + Blin*(y(i)+q(i) ) );
    ym = Clin*x;
    r(i+1) = rho*ym;
end

figure
plot(t, y(2:end))
hold on
plot(t, r(2:end))
% ylim([0, 1])
% xlim([0, 0.1])

a_lin 
b_nonlin

b_lin
a_nonlin