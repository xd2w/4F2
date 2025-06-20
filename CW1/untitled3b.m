close all;
figure
G_nom = tf(1, [1, 1, 1]);
% S = 1/G_nom;
dcp = 0.075;
dcv = 0.1;
w_cr = sqrt(dcp/(dcp+dcv));
s = tf('s');
W1 = makeweight(100, 6*pi, 0.01);

% H = makeweight(0.1, w_cr, 1.1);

D1 = ((s^2 + s+ 1)*(s^2 + (1+dcv)*s +1 + dcp));
D2 = ((s^2 + s+ 1)*(s^2 + (1+dcv)*s +1 + dcp));

W21 = (-dcp*(s^2+s+1) - dcv*s +dcp)/((s^2 + s+ 1)*(s^2 + (1+dcv)*s +1  -dcp));
W22 = (dcp*(s^2+s+1) + dcv*s -dcp)/((s^2 + s+ 1)*(s^2 + (1-dcv)*s +1 + dcp));
W23 = W22 - W21;
W2 = 0.3/(1/W21 + 1/W22);
W2p = (0.6*(s) / ((s +1)^3));
display(W23);

% W2 = sqrt(w^2*(dcp^2*w^2 + (dcp+dcv)^2)/(w^4 - w^2 +1)/(((1-dcv)^2)*w^2 + (1+dcp-w^2)^2) )
    % (((1-heaviside(w-1i*w_cr))/(((1-dcv)^2)*w^2 + (1+dcp-w^2)^2)) + (heaviside(w-1i*w_cr)/(((1+dcv)^2)*w^2 + (1-dcp-w^2)^2)));
W3 = makeweight(0.01, 500*pi, 100);
P = augw(G_nom,W1,[] ,W3);
[K,CL,gamma] = hinfsyn(P, 1, 1);
gamma


[K,CL,gamma] = mixsyn(G_nom, W1, W23, W3);
gamma

figure
% hinfsyn();
% mixsyn();
bodemag(W21, W22, W23)
legend()

figure
bodemag(W1, W2, W2p, W3)
legend()