%% Équations analytiques du modèle
function dx = synchronous_machine(x, Rs, D, J, p, Vs, Td, ws, If, id_interp, iq_interp, Tem_interp)
% variables d'état
% dx(1) = dflux_d 
% dx(2) = dflux_q
% dx(3) = ddelta
% dx(4) = ddw
% dx(5) = dtheta
dx = zeros(5,1);

% fréquence électrique de la machine synchrone (rad/s)
w = ws + x(4); 

% calcul de l'angle du rotor 
dx(5) = w/p;

% calcul des tensions d'axe direct et en quadrature
vd = -Vs*sin(x(3));
vq =  Vs*cos(x(3)); 

% calcul des courants d'axe direct et en quadrature
% et du couple électromagnétique, en normalisant theta dans
% l'intervalle
theta = mod(rad2deg(x(5)),60);
id = id_interp(x(1), x(2), If, theta);
iq = iq_interp(x(1), x(2), If, theta);Tem = Tem_interp(id, iq, If, theta);

% équations de transition d'état
dx(1) = vd - Rs*id + w*x(2);
dx(2) = vq - Rs*iq - w*x(1);
dx(3) = x(4);
dx(4) = (p*(Td - Tem)-D*x(4))/J;

end