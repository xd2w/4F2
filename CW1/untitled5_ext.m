close all;
tau = 0.05;
% X = tf(1, [tau^2, 2*tau, 1]);
s = tf('s');
X = 1/(tau*s + 1)^2;
W = 1/(0.5*s +1)^4;

G = [[X*s^2  -X  -X  0];
     [X (s^(-2))*(1/X - X) -(s^(-2))*X (s^(-2))];
     [W*(1-X) -W*(s^(-2))*(1/X - X) W*(s^(-2))*X W*(-s^(-2))];
     [W*(s^2)*X W*(1-X) -W*X 0]];

% makeStateSpace(G)
Gss = ss(makeStateSpace(G), 'minimal');

% sigma(Gss);

%%%
A = Gss.A;
B = Gss.B;
C = Gss.C;
D = Gss.D;
Br = B(:, 1:2);
Bu = B(:, 3:4);
% Dr = D(:, 1:2);
% Du = D(:, 3:4);
Dry = D(1:2, 1:2);
Duy = D(1:2, 3:4);
Dre = D(3:4, 1:2);
Due = D(3:4, 3:4);
Cy = C(1:2, :);
Ce = C(3:4, :);


[n, n] = size(A);
% gamma = 100;


% cvx_begin sdp
%     variable P(n , n) symmetric
%     % variable gamma
%     minimize(trace(P))
%     P >= 0
%     % gamma >= 0
%     A'*P + P*A <= 0
% cvx_end

n=12;
gamma = 0;

% cvx_begin sdp
%     variable Y(n, n) symmetric
%     variable Z(2, n)
%     variable gamma1
%     minimize(gamma1)
%     LMI1 = Y>0
%     LMI2 = gamma1 > 0.1
%     LMI3 = [[Y*A' +  A*Y + Z'*Bu' + Bu*Z , Br, Y*Cy' + Z'*Duy' ];
%                 [ Br', -gamma1*eye(2) , Dry'];
%                 [ Cy*Y + Duy*Z , Dry , -gamma1*eye(2)]] < 0;
% cvx_end
% 
% K1 = Z/Y

cvx_begin sdp
    variable Y(n, n) symmetric
    variable Z(2, n)
    variable gamma1
    minimize(gamma1)
    LMI1 = Y>0
    LMI2 = gamma1 > 0.1
    LMI3 = [[Y*A' +  A*Y + Z'*Bu' + Bu*Z , Br, Y*Ce' + Z'*Due' ];
                [ Br', -gamma1*eye(2) , Dre'];
                [ Ce*Y + Due*Z , Dre , -gamma1*eye(2)]] < 0;
cvx_end

K = Z/Y;

size(K)
s = tf('s');
% G_feedback = (Ce + Due*K)/(eye(n)*s-(A+Bu*K))*Br
G_feedback = ss((A+Bu*K), Br, (Ce + Due*K), Dre)
K_full = [zeros(2, n); K];
G_full_fb = ss(A+B*K_full, B, C+D*K_full, D);

% Gre = ss(A, Bu, Ce, Due);
% Gfb = feedback()
figure 
% step(G_feedback)
step(G_full_fb)

figure
t = 0:1e-2:10;
u = [sin(10*t); sin(10*t)];
[y_out, time] = lsim(G_feedback, u, t);
plot(time, y_out)

figure
t = 0:1e-5:10;
u = [sin(1e4*t); sin(1e4*t)];
[y_out, time] = lsim(G_feedback, u, t);
plot(time, y_out)

figure
t = 0:1e-5:10;
u = [sin(10*t); sin(1e4*t)];
[y_out, time] = lsim(G_feedback, u, t);
plot(time, y_out)


figure
t = 0:1e-5:10;
u = [sin(1e4*t); sin(10*t)];
[y_out, time] = lsim(G_feedback, u, t);
plot(time, y_out)

figure
sigma(G_full_fb)

% figure
hinfnorm(G_full_fb)

tfs = tf(G_full_fb);
Gq1_fh = tfs(1, 1);
Gfe_q2 = tfs(2, 2);

hinfnorm(Gq1_fh)

hinfnorm(Gfe_q2)

figure 
bodemag(Gq1_fh, Gfe_q2)

% figure
% step(Gss)

