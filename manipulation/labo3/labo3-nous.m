clear all;close all;clc
h = 0.1;
nH = 1;
dH = conv([1 0], [1,0.5]);
Hc = tf(nH, dH);
Hd = c2d(Hc, h, 'zoh');

% R(z) = [K(Z-p)(Z-n)]/[(Z-z)(Z-1)]
[z,p,~] = zpkdata(Hd, 'v');
p = p(2);

[pc] = dtr2ord2o(1, 5, 5); % trouve les poles pd1 et pd2
pd = exp(pc*h);            % discretise les poles
K = 2-pd(1)-pd(2);         % calcul le gain d'Evans unite
n = (1-pd(1)*pd(2))/K;     % calcul du zero additionnel

Z = tf('z', h);
F = K*(Z-n)/((Z-1)^2+K*(Z-n));    % boucle fermee
R = K*(Z-p)*(Z-n)/((Z-z)*(Z-1));  % regulateur

% autre maniere de fermer la boucle
Bo = R*Hd;
Bf = feedback(Bo,1);
Bf = minreal(Bf);
