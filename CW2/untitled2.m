s = tf('s');

A = [
    -10,    0,  10;
      0, -100, 100;
      0,    0,  -1;
    ];

B = [0; 0; 1];
C = [0, 0, 1];
D= 0;
[n, ~] = size(A);

lambda = 20;
ep = 1e-50;
gamma = 2
cvx_begin sdp
    cvx_solver SDPT3

    variable Y(n, n) symmetric
    variable Z(1, n)
    % variable gamma
    % variable ep
    % minimize(ep)
    % LMI1 = Y>0
    % LMI2 = gamma < 2
    LMI3 = ep > 0
    % LMI4 = Y*A' +  A*Y + Z'*B' + B*Z + 2*lambda*Y + ep*eye(n) < 0
    % LMI5 = Y*A' +  A*Y + 2*lambda*Y + ep*eye(n) < 0
    LMI6 = Y*A' +  A*Y + Z'*B' + B*Z + ep*eye(n) < 0
    LMI4 = [[Y*A' +  A*Y + Z'*B' + B*Z + 2*lambda*Y + ep*eye(n), B, Y*C'];
                [ B', -gamma , 0];
                [ C*Y , 0 , -gamma]] < 0;
    LMI5 = [[Y*A' +  A*Y + 2*lambda*Y + ep*eye(n), B, Y*C'];
                [ B', -gamma , 0];
                [ C*Y , 0 , -gamma]] < 0;
cvx_end

K = Z/Y;
% Y

nt = 10;
N = 10000;
t = linspace(0, nt, N);
dt = nt/(N);

r = ones(1, N)*0.5; 
% r = sin(t);
y = zeros(1, N);

z = [0; 0; 0];

for i=1:N
    z = z + dt* ((A*z) + B*(tanh(K*z) + r(i)));
    y(i) = z(3);
end

figure
plot(t, y)
hold on;

E = sum(y.*y)/sum(r.*r)

z = [1; 1; 1];

for i=1:N
    z = z + dt* ((A*z) + B*(tanh(K*z) + r(i)));
    y(i) = z(3);
end

plot(t, y)
xlabel("Time")
ylabel("\delta y")
saveas(gcf, './figs/1iii_r=0-5', 'epsc')

% E = sum(y.*y)/sum(r.*r)
% plot(t(1:end),E);


r = ones(1, N)*1; 
% r = sin(t);
y = zeros(1, N);

z = [0; 0; 0.001];

for i=1:N
    z = z + dt* ((A*z) + B*(tanh(K*z) + r(i)));
    y(i) = z(3);
end

figure
plot(t, y)
hold on;

E = sum(y.*y)/sum(r.*r)

z = [1; 1; 1];

for i=1:N
    z = z + dt* ((A*z) + B*(tanh(K*z) + r(i)));
    y(i) = z(3);
end

plot(t, y)
xlabel("Time")
ylabel("\delta y")
saveas(gcf, './figs/1iii_r=1', 'epsc')

% E = sum(y.*y)/sum(r.*r)
% plot(t(1:end),E);


eig(Y)

% eig(A)