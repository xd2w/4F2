close all;
s = tf('s');
tau = 0.05;
X = 1/(tau*s + 1)^2;

G = [[1/s, X/s 1/s, 0];
    [X/s, 1/s, 0, 1/s];
    [(1-X)/(s^2), (X-1)/(s^2), s^(-2), -s^(-2)];
    [(1-X)/s, (X-1)/s, 1/s, -1/s]];

Gss  = ss(makeStateSpace(G), 'minimal');
isPassive(Gss)
step(Gss)

% sigma(Gss)
% hinfnorm(Gss)

% isproper(minreal(Gss))

% figure
% [Gss, Gred_info] = ncfmr(minreal(Gss), 5)
% reducespec(Gss, 'ncf')

A = Gss.A
B = Gss.B;
C = Gss.C;
D = Gss.D;
Bf = B(:, 1:2);
Bu = B(:, 3:4);
% Df = D(:, 1:2);
% Du = D(:, 3:4);
Dfq = D(1:2, 1:2);
Duq = D(1:2, 3:4);
Dfe = D(3:4, 1:2);
Due = D(3:4, 3:4);
Cq = C(1:2, :);
Ce = C(3:4, :);


[n, n] = size(A);

ep = 0;
gamma =0;

cvx_begin sdp
    variable Q(n, n) symmetric
    variable Y(n, n) symmetric
    variable Y(2, n)
    variable Z(2, n)
    variable ep
    variable gamma
    minimize(ep)
    LMI1 = Q>0
    LMI2 = ep > 0.
    LMI3 = [Q*A' + A*Q + Bu*Y + Y'*Bu', Bf - Q*Cq' - Y'*Duq';
            Bf' - Cq*Q - Duq*Y, -2*ep*eye(2) - (Dfq + Dfq')] <= 0;

    LMI4 = Yp>0
    LMI5 = gamma > 0.1
    LMI6 = [[Y*A' +  A*Y + Z'*Bu' + Bu*Z , Br, Y*Ce' + Z'*Due' ];
                [ Br', -gamma*eye(2) , Dre'];
                [ Ce*Y + Due*Z , Dre , -gamma*eye(2)]] < 0;
cvx_end

K = Y/Q;

K_full = [zeros(2, n); K];

G_feedback = ss((A+Bu*K), Bf, (Cq + Duq*K), Dfq);
Gfb_full = ss((A+B*K_full), B, (C + D*K_full), D);
% Gfe = ss((A+Bf*K), Bf, (Ce + Due*K), Due);
isPassive(G_feedback)

tfs = tf(Gfb_full);

hinfnorm(tfs(3:4, 1:2))

figure
sigma(Gfb_full)

figure
% sigma(tfs(1:2, 3:4))
sigma(tfs(3:4, 1:2))
bode(tfs(1:2, 1:2))
hinfnorm(tfs(3:4, 1:2))

figure
step(G_feedback)

figure
step(Gfb_full)