function [] = Init()
%% Constantes de la machine synchrone

Vs = 70;
fs = 50;
ws = 2*pi*fs;
Rs = 0.05;
Xq = 0.7874;
Xaf = 3.6742;
Xd = 1.5903;
D = 10;
J = 0.5;
p = 2;
fluxd0 = 0.2183;
fluxq0 = 0.0525;
delta0 = 0.2250;
dw0 = 0;
theta0 = 0;
%% Données de l'analyse FEM de la machine

luts = load("LUT_FINAL_SIMULINKv2.mat"); % données de courant et de couple
grad = load("GRADIENTS.mat");            % dérivées numériques 
Sac = load("S.mat");                     % points de fonctionnement de la machine
Pac = Sac.map_P;
Qac = Sac.map_Q;

%% Création des interpolateurs des cartographies de courant et de couple
id_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, luts.Id, 'linear', 'linear'); % interpolateur calcul du courant direct
iq_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, luts.Iq, 'linear', 'linear'); % interpolateur calcul du courant en quadrature
Tem_interp = griddedInterpolant({luts.vId, luts.vIq, luts.vIf, luts.vTheta}, luts.Tem, 'linear', 'linear');   % interpolateur calcul du couple

% interpolateur des gradients 
did_fd_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, grad.Gidfd, 'linear', 'linear');       % id/fluxd
did_fq_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, grad.Gidfq, 'linear', 'linear');       % id/fluxq
did_theta_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, grad.Gidtheta, 'linear', 'linear'); % id/theta
diq_fd_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, grad.Giqfd, 'linear', 'linear');       % iq/fluxd
diq_fq_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, grad.Giqfq, 'linear', 'linear');       % iq/fluxq
diq_theta_interp = griddedInterpolant({luts.vPsid, luts.vPsiq, luts.vIf, luts.vTheta}, grad.Giqtheta, 'linear', 'linear'); % iq/theta
dTem_id_interp = griddedInterpolant({luts.vId, luts.vIq, luts.vIf, luts.vTheta}, grad.GTemid, 'linear', 'linear');         % Tem/fluxd
dTem_iq_interp = griddedInterpolant({luts.vId, luts.vIq, luts.vIf, luts.vTheta}, grad.GTemiq, 'linear', 'linear');         % Tem/fluxq
dTem_theta_interp = griddedInterpolant({luts.vId, luts.vIq, luts.vIf, luts.vTheta}, grad.GTemtheta, 'linear', 'linear');   % Tem/theta

variaveis = whos; 
for i = 1:length(variaveis)
    assignin('base', variaveis(i).name, eval(variaveis(i).name));
end

end