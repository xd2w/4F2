cla;
delcp_range = -0.075:0.0075:	0.075;
delcv_range = -0.1:0.01:0.1;

w = logspace(-2, 2, 100);
gain = zeros(size(w));
% % s = 1i * w;

% s = tf('s');

figure
hold on;
axes('XScale', 'log', 'YScale', 'log')
for delcp = delcp_range
    for delcv = delcv_range
        new_gain = ( ...
            ((delcp^2).*(w.^4) + (w.^2).*(delcp-delcv)^2) ...
            ./ ...
            ((1-w.^2+w.^4).*(((1+delcv)^2)*(w.^2) + (1+delcp-w.^2).^2)) ...
            ).^0.5;
        % loglog(w, new_gain);
        gain = max(gain, 10*new_gain.^2);
        hold on;
    end
end
legend
plot(w, gain);
hold on;

n = 21;

delcp = delcp_range(n);
delcv = delcv_range(1);
gain = ( ...
            ((delcp^2).*(w.^4) + (w.^2).*(delcp-delcv)^2) ...
            ./ ...
            ((1-w.^2+w.^4).*(((1+delcv)^2)*(w.^2) + (1+delcp-w.^2).^2)) ...
            ).^0.5;
plot(w, 10*gain.^2);


delcp = delcp_range(1);
delcv = delcv_range(n);
gain = ( ...
            ((delcp^2).*(w.^4) + (w.^2).*(delcp-delcv)^2) ...
            ./ ...
            ((1-w.^2+w.^4).*(((1+delcv)^2)*(w.^2) + (1+delcp-w.^2).^2)) ...
            ).^0.5;
plot(w, 10*gain.^2);

dcp = 0.075;
dcv = 0.1;

D1 = ((s^2 + s+ 1)*(s^2 + (1+dcv)*s +1 + dcp));
D2 = ((s^2 + s+ 1)*(s^2 + (1+dcv)*s +1 + dcp));

W21 = (-dcp*(s^2+s+1) - dcv*s +dcp)/((s^2 + s+ 1)*(s^2 + (1+dcv)*s +1  -dcp));
W22 = (dcp*(s^2+s+1) + dcv*s -dcp)/((s^2 + s+ 1)*(s^2 + (1-dcv)*s +1 + dcp));

W2 = 1/(1/W21 + 1/W22);
W2p = (1.6*(s) / (s +0.9)^2);
figure
bode(W2)


% hod on
% hold off
% 
% G = 1 + tf(1, [1 , 1, 1]);
% bode(G);

