close all;
m=1;
cv=1+0*0.075;
cp=1+0*0.1;
G_nom = tf(cp, [m, cv, cp]);

figure(5);
% bode(G_nom);
margin(G_nom)
% disp('inf gain margin  90 phase margin');

figure(3);
nyquist(G_nom);

figure(2);
rlocus(G_nom);


figure(10)
t = 0:1e-4:10e-1;
omega = 2*pi*1e2;
u = sin(omega*t);
% d = tf(omega, [1, 0, omega^2]);

% step = tf(1, [1, 0]);
for k=1:1:9
    % k= -k;
    % y = -((k*G_nom)/(1+k*G_nom))*d;
    % step(((k*G_nom)/(1+k*G_nom)))
    [yout, time] = lsim(((k*G_nom)/(1+k*G_nom)), u, t);
    plot(time, yout, 'DisplayName', num2str(k))
    hold on;
end

cp = ureal('cp', 1, 'PlusMinus', 0.075);
cv = ureal('cv', 1, 'PlusMinus', 0.1);
s = tf('s');
G_nom_u = cp/(s^2 + cv*s + cp);

figure
nyquist(usample(G_nom_u, 100))