clc; clear; close all;

%% Iniatilization
Init();

%% Simulation en régime permanent 

load("Po.mat") % carregando diferentes pontos de operação, Td_If (7x2) : [Td, If]
% Po(1,1) = 85.3373; Po(1,2) = 650.38;      Qac > 0 
% Po(2,1) = 1622.16; Po(2,2) = 460.719;     Qac > 0
% Po(3,1) = 2830.62; Po(3,2) = -823.729;    Qac < 0
% Po(4,1) = 3386.87; Po(4,2) = -3075.32;    Qac < 0
% Po(5,1) = 2690.85; Po(5,2) = -5152.82;    Qac < 0 
% Po(6,1) = 1781.44; Po(6,2) = -2370.5;     Qac < 0
% Po(7,1) = 294.574; Po(7,2) = -1859.97;    Qac < 0
% Po(8,1) = 1719; Po(8,2) = -1001.15;    % Qac < 0
point = 8; % ponto de operação escolhido

dt = 1e-4; % pas de simulation
tf = 2; % temps de simulation
tav = 0:dt:tf/2; % temps avant la pertubation (changement de Td)
tap = tav(end)+dt:dt:tf; % temps après la pertubation (changement de Td)
% tav = 0:dt:tf/2; % temps avant la pertubation (court circuit)
% tpd = tav(end)+dt:dt:tf/2+0.01; % temps après la pertubation (court circuit)
% tap = tpd(end)+dt:dt:tf; % temps avant la pertubation (court circuit)

If = Td_If(point,2);
Td = Td_If(point,1);

x0 = [fluxd0, fluxq0, delta0, dw0, theta0]'; % point de départ des états
[t1, x1] = ode15s(@(t,x)synchronous_machine(x, Rs, D, J, p, Vs, Td, ws, If, id_interp, iq_interp, Tem_interp), tav, x0); 
[t2, x2] = ode15s(@(t,x)synchronous_machine(x, Rs, D, J, p, Vs, 0.5*Td, ws, If, id_interp, iq_interp, Tem_interp), tap, x1(end,:)); 
% [t1, x1] = ode15s(@(t,x)synchronous_machine(x, Rs, D, J, p, Vs, Td, ws, If, id_interp, iq_interp, Tem_interp), tav, x0); 
% [t2, x2] = ode15s(@(t,x)synchronous_machine(x, Rs, D, J, p, 0.8*Vs, Td, ws, If, id_interp, iq_interp, Tem_interp), tpd, x1(end,:)); 
% [t3, x3] = ode15s(@(t,x)synchronous_machine(x, Rs, D, J, p, Vs, Td, ws, If, id_interp, iq_interp, Tem_interp), tap, x2(end,:)); 

x = [x1; x2]; % changement de Td
t = [t1; t2];
% x = [x1; x2; x3]; % court circuit
% t = [t1; t2; t3];

vd = -Vs.*sin(x(:,3));
vq = Vs.*cos(x(:,3));
id = zeros(1,length(t));
iq = zeros(1,length(t));

for i=1:length(t)
    id(1,i) = id_interp(x(i,1), x(i,2), If, mod(rad2deg(x(i,5)),60));
    iq(1,i) = iq_interp(x(i,1), x(i,2), If, mod(rad2deg(x(i,5)),60));
end

save("id_meas.mat", "id");
save("iq_meas.mat", "iq");
save("vd_meas.mat", "vd");
save("vq_meas.mat", "vq");
%% Tracé en régime permanent
% tracé des états
figure;
plot(t,x(:,1), 'LineWidth', 1.5);
title('Flux direct');
xlabel('Temps (sec)');
ylabel('Flux (Wb)');

figure;
plot(t,x(:,2), 'LineWidth', 1.5);
title('Flux en quadrature');
xlabel('Temps (sec)');
ylabel('Flux (Wb)');

figure;
plot(t,rad2deg(x(:,3)), 'LineWidth', 1.5); 
title('Angle interne');
xlabel('Temps (sec)');
ylabel('\delta (Deg)');

figure;
plot(t,x(:,4), 'LineWidth', 1.5); 
title('Variation de vitesse');
xlabel('Temps (sec)');
ylabel('\Delta\omega (rad/s)');

figure;
plot(t,mod(rad2deg(x(:,5)),60), 'LineWidth', 1.5); 
title('Position du rotor');
xlabel('Temps (sec)');
ylabel('\theta (Deg)');

% tracé des puissances
vd = -Vs.*sin(x(:,3));
vq = Vs.*cos(x(:,3));
id = zeros(1,length(t));
iq = zeros(1,length(t));

for i=1:length(t)
    id(1,i) = id_interp(x(i,1), x(i,2), If, mod(rad2deg(x(i,5)),60));
    iq(1,i) = iq_interp(x(i,1), x(i,2), If, mod(rad2deg(x(i,5)),60));
end
[P, Q] = S(vd, vq, id, iq);

figure;
subplot(2,1,1);
plot(t, P, 'LineWidth', 1.5); title('Puissance Active'); xlabel('Temps (sec)'); ylabel('P (W)'); 
legend('Puissance Actif');
subplot(2,1,2);
plot(t, Q, 'LineWidth', 1.5); title('Puissance Réactive'); xlabel('Temps (sec)'); ylabel('Q (VAr)'); 
legend('Puissance Réactive');

%% Analyse de la puissance active en fonction du couple mécanique
% Et puissance réactive en fonction du courant d'excitation
Td_f = 100;
If_f = 40;

Td_sim = 0:Td_f/50:Td_f; % variation du coupale mécanique
If_sim = 5:(If_f-5)/50:If_f; % variation du courant d'excitation

dt = 1e-4; % pas de simulation
tf = 0.5; % temps de simulation
t = 0:dt:tf; % temps de simulation

x0 = [fluxd0, fluxq0, delta0, dw0, theta0]'; % point de départ des états

map_P = zeros(length(Td_sim), length(If_sim));
map_Q = zeros(length(Td_sim), length(If_sim));
%% Génération des cartographies de puissance en fonction de Td et If

for i=1:length(Td_sim)
    for j=1:length(If_sim)
        fprintf("%d, %d\n", i, j);
        [t, x] = ode15s(@(t,x)synchronous_machine(x, Rs, D, J, p, Vs, Td_sim(i), ws, If_sim(j), id_interp, iq_interp, Tem_interp), t, x0); 
        vd = -Vs.*sin(x(:,3));
        vq = Vs.*cos(x(:,3));
        id = zeros(1,length(t));
        iq = zeros(1,length(t));
        for k=1:length(t)
            id(1,k) = id_interp(x(k,1), x(k,2), If_sim(j), mod(rad2deg(x(k,5)),60));
            iq(1,k) = iq_interp(x(k,1), x(k,2), If_sim(j), mod(rad2deg(x(k,5)),60));
        end
        [P, Q] = S(vd, vq, id, iq);
        map_P(i,j) = mean(P);
        map_Q(i,j) = mean(Q);
    end
end

save("S.mat", "map_P", "map_Q");

%% Courbes de la Puissance Active par rapport à la Puissance Réactive (points de fonctionnement de la machine)
figure;
surf(Td_sim, If_sim, Pac);

figure;
surf(Td_sim, If_sim, Qac);

figure;
scatter(reshape(Qac,[],1), reshape(Pac,[],1));
title('Courbe P x Q');
xlabel('Q [VAr]');
ylabel('P [W]');

%% Détermination de deux points de fonctionnement pour les cas où Pac > 0 et Qac > 0 ou Qac < 0
Po = zeros(8,2); % ponto de operação (P, Q)

Po(1,1) = 85.3373; Po(1,2) = 650.38;     % Qac > 0
Po(2,1) = 1622.16; Po(2,2) = 460.719;    % Qac > 0
Po(3,1) = 2830.62; Po(3,2) = -823.729;    % Qac < 0
Po(4,1) = 3386.87; Po(4,2) = -3075.32;     % Qac < 0
Po(5,1) = 2690.85; Po(5,2) = -5152.82;    % Qac < 0
Po(6,1) = 1781.44; Po(6,2) = -2370.5;     % Qac < 0
Po(7,1) = 294.574; Po(7,2) = -1859.97;    % Qac < 0
Po(8,1) = 1719; Po(8,2) = -1001.15;    % Qac < 0

Td_If = zeros(8,2); % pares Td If para cada ponto de operação (Td, If)

for i=1:8
    erroP = abs(Pac - Po(i,1));
    erroQ = abs(Qac - Po(i,2));
    [~, idx1P] = min(erroP(:));     % index P1 
    [~, idx1Q] = min(erroQ(:));     % index Q1
    if idx1P==idx1Q
        [idx1r, idx1c] = ind2sub(size(Pac), idx1P);
    end

    Td_If(i,1) = Td_sim(idx1r); % Td du point de fonctionnement choisi
    Td_If(i,2) = If_sim(idx1c); % If du point de fonctionnement choisi
end

save("Po.mat", "Td_If");

%% Filtres
% chargement des états réels de la machine pour calculer les sorties réelles
id = load("id_meas.mat").id;
iq = load("iq_meas.mat").iq;
vq = load("vq_meas.mat").vq';

n = length(id);

rho = 1;
Ke = 1;
iex = Xaf*If/Ke;

% Valores reais das saídas (sem ruído)
Pac_meas = (iq.*(vq - Xq*id));
Qac_meas = (-Xq*iq.^2 - id.*vq);
Is_2_meas = (id.^2 + iq.^2);
Vs_2_meas = (Xq^2*iq.^2 + vq.^2);
y5 = tanh(rho./iex.*(vq -Xd*id));

% Inserindo ruído Gaussiano (média zero) com intensidade B
bu = 0.05;
Pac_meas_ = Pac_meas + bu*std(Pac_meas)*randn(size(id));
Qac_meas_ = Qac_meas + bu*std(Qac_meas)*randn(size(id));
Is_2_meas_ = Is_2_meas + bu*std(Is_2_meas)*randn(size(id));
Vs_2_meas_ = Vs_2_meas + bu*std(Vs_2_meas)*randn(size(id));
y5_meas = tanh(rho./iex.*(vq - Xd*id)) + bu*std(y5)*randn(size(id));

meas = [Pac_meas_; Qac_meas_; Is_2_meas_; Vs_2_meas_; y5_meas];

x0 = [6, 10, 50];
P0 = eye(3);
Q = diag([1e-2, 1e-2, 1e-5]);
R = diag([
    var(Pac_meas_ - Pac_meas), ...
    var(Qac_meas_ - Qac_meas), ...
    var(Is_2_meas_ - Is_2_meas), ...
    var(Vs_2_meas_ - Vs_2_meas), ...
    1e-4
]);

x = x0';
P = P0;
%% Utilisation de l'EKF pour une simulation

[x_est_EKF, y_est_EKF] = EKF(Q, R, x, P, rho, Ke, Xd, Xq, Xaf, If, meas, n);

%% comparaison des états réels et ceux estimés par l'EKF
% calcul du RMSE
RMSE = zeros(3,n);
for i=1:n
    RMSE(:,i) = rmse(x_est_EKF(:,1:i),[id(1:i); iq(1:i); vq(1:i)], 2);
end

%%
figure;
subplot(2,1,1);
plot(t, id, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, x_est_EKF(1,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Courant d''axe direct I_d');
xlabel('Temps (sec)');
ylabel('I_d');
subplot(2,1,2)
plot(t, RMSE(1,:), 'LineWidth', 1.5, 'Color', 'red');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, iq, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, x_est_EKF(2,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Courant d''axe en quadrature I_q');
xlabel('Temps (sec)');
ylabel('I_q');
subplot(2,1,2)
plot(t, RMSE(2,:), 'LineWidth', 1.5, 'Color', 'red');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, vq, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, x_est_EKF(3,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Tension d''axe direct');
xlabel('Temps (sec)');
ylabel('v_d');
subplot(2,1,2)
plot(t, RMSE(3,:), 'LineWidth', 1.5, 'Color', 'red');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, Pac_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_EKF(1,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Puissance active');
xlabel('Temps (sec)');
ylabel('P_{ac}');
subplot(2,1,2)
plot(t, Qac_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_EKF(2,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Puissance réactive');
xlabel('Temps (sec)');
ylabel('Q_{ac}');

figure;
subplot(3,1,1);
plot(t, Is_2_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_EKF(3,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Courant terminal');
xlabel('Temps (sec)');
ylabel('I_s^2');
subplot(3,1,2)
plot(t, Vs_2_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_EKF(4,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Tension terminale');
xlabel('Temps (sec)');
ylabel('V_s^2');
% subplot(3,1,3)
% plot(t, y5, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
% plot(t, y_est_EKF(5,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
% title('Cinquième sortie');
% xlabel('Temps (sec)');
% ylabel('...');
%% Utilisation de l'UKF

[x_est_UKF, y_est_UKF] = UKF(x, P, Q, R, Xd, Xq, Xaf, rho, Ke, If, meas, n);

%% comparaison des états réels et ceux estimés par l'UKF
% calcul du RMSE
RMSE = zeros(3,n);
for i=1:n
    RMSE(:,i) = rmse(x_est_UKF(:,1:i),[id(1:i); iq(1:i); vq(1:i)], 2);
end

%%
figure;
subplot(2,1,1);
plot(t, id, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, x_est_UKF(1,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Courant d''axe direct I_d');
xlabel('Temps (sec)');
ylabel('I_d');
subplot(2,1,2)
plot(t, RMSE(1,:), 'LineWidth', 1.5, 'Color', 'red');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, iq, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, x_est_UKF(2,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Courant d''axe en quadrature I_q');
xlabel('Temps (sec)');
ylabel('I_q');
subplot(2,1,2)
plot(t, RMSE(2,:), 'LineWidth', 1.5, 'Color', 'red');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, vq, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, x_est_UKF(3,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Tension d''axe direct');
xlabel('Temps (sec)');
ylabel('v_d');
subplot(2,1,2)
plot(t, RMSE(3,:), 'LineWidth', 1.5, 'Color', 'red');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, Pac_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_UKF(1,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Puissance active');
xlabel('Temps (sec)');
ylabel('P_{ac}');
subplot(2,1,2)
plot(t, Qac_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_UKF(2,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Puissance réactive');
xlabel('Temps (sec)');
ylabel('Q_{ac}');

figure;
subplot(3,1,1);
plot(t, Is_2_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_UKF(3,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Courant terminal');
xlabel('Temps (sec)');
ylabel('I_s^2');
subplot(3,1,2)
plot(t, Vs_2_meas, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_UKF(4,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Tension terminale');
xlabel('Temps (sec)');
ylabel('V_s^2');
subplot(3,1,3)
plot(t, y5, 'LineWidth', 1.5, 'Color', 'blue'); hold on;
plot(t, y_est_UKF(5,:), 'LineWidth',  0.8, 'Color', 'black', 'LineStyle', '--');
title('Cinquième sortie');
xlabel('Temps (sec)');
ylabel('...');

%% Expériences de Monte Carlo pour l'EKF et l'UKF
% création des valeurs initiales 
X = 100; % nombre d'expériences
m = 1.5;  % indice de variation (0 : 0% | 1 : 100%)

% vecteur de valeurs initiales aléatoires
id_0 = id(1) + m*id(1)*(rand(1,X) - 0.5)*2;
iq_0 = iq(1) + m*iq(1)*(rand(1,X) - 0.5)*2;
vq_0 = vq(1) + m*vq(1)*(rand(1,X) - 0.5)*2;

% création des vecteurs de résultats
x_EKF_MC = zeros(3, n, X);
y_EKF_MC = zeros(4, n, X);
% x_UKF_MC = zeros(3, n, X);
% y_UKF_MC = zeros(5, n, X);

for i=1:X
    % attribution des différentes valeurs initiales pour l'EKF
    x = [id_0(i), iq_0(i), vq_0(i)]';
    P = P0;
    [x_EKF_MC(:,:,i), y_EKF_MC(:,:,i)] = EKF(Q, R, x, P, rho, Ke, Xd, Xq, Xaf, If, meas, n);

    % réinitialisation de la matrice et attribution des valeurs pour l'UKF
    % x = [id_0(i), iq_0(i), vq_0(i)]';
    % P = P0;
    % [x_UKF_MC(:,:,i), y_UKF_MC(:,:,i)] = UKF(x, P, Q, R, Xd, Xq, Xaf, rho, Ke, If, meas, n);

    fprintf('Progress %.2f%%\n', i/X*100);
 
end

%% Comparaison entre les filtres

% états
% calcul du RMSE
RMSE = zeros(3,n,2);
for i=1:n
    RMSE(:,i,1) = rmse(x_est_EKF(:,1:i),[id(1:i); iq(1:i); vq(1:i)], 2);
    RMSE(:,i,2) = rmse(x_est_UKF(:,1:i),[id(1:i); iq(1:i); vq(1:i)], 2);
end

% tracé
figure;
subplot(2,1,1);
plot(t, id, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, x_est_EKF(1,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, x_est_UKF(1,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
legend('Réel', 'EKF', 'UKF');
title('Courant d''axe direct I_d');
xlabel('Temps (sec)');
ylabel('I_d');
subplot(2,1,2)
plot(t, RMSE(1,:,1), 'LineWidth', 1.5, 'Color', 'red'); hold on;
plot(t, RMSE(1,:,2), 'LineWidth', 1.5, 'Color', 'magenta');
legend('EKF', 'UKF');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, iq, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, x_est_EKF(2,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, x_est_UKF(2,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
legend('Réel', 'EKF', 'UKF');
title('Courant d''axe en quadrature I_q');
xlabel('Temps (sec)');
ylabel('I_q');
subplot(2,1,2)
plot(t, RMSE(2,:,1), 'LineWidth', 1.5, 'Color', 'red'); hold on;
plot(t, RMSE(2,:,2), 'LineWidth', 1.5, 'Color', 'magenta');
legend('EKF', 'UKF');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(2,1,1);
plot(t, vq, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, x_est_EKF(3,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, x_est_UKF(3,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
legend('Réel', 'EKF', 'UKF');
title('Tension d''axe en quadrature');
xlabel('Temps (sec)');
ylabel('v_d');
subplot(2,1,2)
plot(t, RMSE(3,:,1), 'LineWidth', 1.5, 'Color', 'red'); hold on;
plot(t, RMSE(3,:,2), 'LineWidth', 1.5, 'Color', 'magenta');
legend('EKF', 'UKF');
title('Erreur');
xlabel('Temps (sec)');
ylabel('RMSE');

% sorties
% calcul du RMSE
RMSEs = zeros(5,n,2);
for i=1:n
    RMSEs(:,i,1) = rmse(y_est_EKF(:,1:i),[Pac_meas(1:i); Qac_meas(1:i); Is_2_meas(1:i); Vs_2_meas(1:i); y5(1:i)], 2);
    RMSEs(:,i,2) = rmse(y_est_UKF(:,1:i),[Pac_meas(1:i); Qac_meas(1:i); Is_2_meas(1:i); Vs_2_meas(1:i); y5(1:i)], 2);
end

% tracé
figure;
subplot(2,1,1);
yyaxis("left")
plot(t, Pac_meas, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, y_est_EKF(1,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, y_est_UKF(1,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
legend('Réel', 'EKF', 'UKF')
title('Puissance active');
xlabel('Temps (sec)');
ylabel('P_{ac}');
yyaxis("right")
plot(t, RMSEs(1,:,1), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '-'); hold on;
plot(t, RMSEs(1,:,2), 'LineWidth', 1.5, 'Color', 'magenta', 'LineStyle', '-');
legend('Réel', 'EKF', 'UKF', 'RMSE EKF', 'RMSE UKF');
xlabel('Temps (sec)');
ylabel('RMSE');
subplot(2,1,2)
yyaxis("left")
plot(t, Qac_meas, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, y_est_EKF(2,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, y_est_UKF(2,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
title('Puissance réactive');
xlabel('Temps (sec)');
ylabel('Q_{ac}');
yyaxis("right")
plot(t, RMSEs(2,:,1), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '-'); hold on;
plot(t, RMSEs(2,:,2), 'LineWidth', 1.5, 'Color', 'magenta', 'LineStyle', '-');
legend('Réel', 'EKF', 'UKF', 'RMSE EKF', 'RMSE UKF');
xlabel('Temps (sec)');
ylabel('RMSE');

figure;
subplot(3,1,1);
yyaxis("left")
plot(t, Is_2_meas, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, y_est_EKF(3,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, y_est_UKF(3,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
title('Courant terminal');
xlabel('Temps (sec)');
ylabel('I_s^2');
yyaxis("right")
plot(t, RMSEs(3,:,1), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '-'); hold on;
plot(t, RMSEs(3,:,2), 'LineWidth', 1.5, 'Color', 'magenta', 'LineStyle', '-');
legend('Réel', 'EKF', 'UKF', 'RMSE EKF', 'RMSE UKF');
xlabel('Temps (sec)');
ylabel('RMSE');
subplot(3,1,2)
yyaxis("left")
plot(t, Vs_2_meas, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, y_est_EKF(4,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, y_est_UKF(4,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
title('Tension terminale');
xlabel('Temps (sec)');
ylabel('V_s^2');
yyaxis("right")
plot(t, RMSEs(4,:,1), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '-'); hold on;
plot(t, RMSEs(4,:,2), 'LineWidth', 1.5, 'Color', 'magenta', 'LineStyle', '-');
legend('Réel', 'EKF', 'UKF', 'RMSE EKF', 'RMSE UKF');
xlabel('Temps (sec)');
ylabel('RMSE');
subplot(3,1,3)
yyaxis("left")
plot(t, y5, 'LineWidth', 1.5, 'Color', 'black'); hold on;
plot(t, y_est_EKF(5,:), 'LineWidth',  0.8, 'Color', 'blue', 'LineStyle', '--'); hold on;
plot(t, y_est_UKF(5,:), 'LineWidth',  0.8, 'Color', 'red', 'LineStyle', '--');
title('Cinquième sortie');
xlabel('Temps (sec)');
ylabel('...');
yyaxis("right")
plot(t, RMSEs(5,:,1), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '-'); hold on;
plot(t, RMSEs(5,:,2), 'LineWidth', 1.5, 'Color', 'magenta', 'LineStyle', '-');
legend('Réel', 'EKF', 'UKF', 'RMSE EKF', 'RMSE UKF');
xlabel('Temps (sec)');
ylabel('RMSE');

%% Graphiques des expériences de Monte Carlo

figure;
for i=1:X
    plot(t, x_EKF_MC(1,:,i));hold on;
end
plot(t, id, 'LineWidth', 2, 'Color', 'black'); % état réel I_d
title('Courant d''axe direct - EKF')
xlabel('Temps (sec)')
ylabel('I_d')

figure;
for i=1:X
    plot(t, x_EKF_MC(2,:,i));hold on;
end
plot(t, iq, 'LineWidth', 2, 'Color', 'black'); % état réel I_q
title('Courant d''axe en quadrature - EKF')
xlabel('Temps (sec)')
ylabel('I_q') 

figure;
for i=1:X
    plot(t, x_EKF_MC(3,:,i));hold on;
end
plot(t, vq, 'LineWidth', 2, 'Color', 'black'); % état réel v_q
title('Tension d''axe en quadrature - EKF')
xlabel('Temps (sec)')
ylabel('I_d')

% figure;
% for i=1:X
%     plot(t, x_UKF_MC(1,:,i));hold on;
% end
% title('Courant d''axe direct - UKF')
% xlabel('Temps (sec)')
% ylabel('I_d')
% 
% figure;
% for i=1:X
%     plot(t, x_UKF_MC(2,:,i));hold on;
% end
% title('Courant d''axe en quadrature - UKF')
% xlabel('Temps (sec)')
% ylabel('I_d')
% 
% figure;
% for i=1:X
%     plot(t, x_UKF_MC(3,:,i));hold on;
% end
% title('Tension d''axe en quadrature - UKF')
% xlabel('Temps (sec)')
% ylabel('I_d')

figure;
for i=1:X
    plot(t, y_EKF_MC(1,:,i));hold on;
end
plot(t, Pac_meas, 'LineWidth', 2, 'Color', 'black'); % sortie réelle P_ac
title('Puissance active - EKF')
xlabel('Temps (sec)')
ylabel('P_{ac}')

figure;
for i=1:X
    plot(t, y_EKF_MC(2,:,i));hold on;
end
plot(t, Qac_meas, 'LineWidth', 2, 'Color', 'black'); % sortie réelle Q_ac
title('Puissance réactive - EKF')
xlabel('Temps (sec)')
ylabel('Q_{ac}')

figure;
for i=1:X
    plot(t, y_EKF_MC(3,:,i));hold on;
end
plot(t, Is_2_meas, 'LineWidth', 2, 'Color', 'black'); % sortie réelle I_s^2
title('Courant terminal - EKF')
xlabel('Temps (sec)')
ylabel('I_s^2')

figure;
for i=1:X
    plot(t, y_EKF_MC(4,:,i));hold on;
end
plot(t, Vs_2_meas, 'LineWidth', 2, 'Color', 'black'); % sortie réelle V_s^2
title('Tension terminale - EKF')
xlabel('Temps (sec)')
ylabel('V_s^2')

% figure;
% for i=1:X
%     plot(t, y_EKF_MC(5,:,i));hold on;
% end
% title('Cinquième sortie - EKF')
% xlabel('Temps (sec)')
% ylabel('...')

% figure;
% for i=1:X
%     plot(t, y_UKF_MC(1,:,i));hold on;
% end
% title('Puissance active - UKF')
% xlabel('Temps (sec)')
% ylabel('P_{ac}')
% 
% figure;
% for i=1:X
%     plot(t, y_UKF_MC(2,:,i));hold on;
% end
% title('Puissance réactive - UKF')
% xlabel('Temps (sec)')
% ylabel('Q_{ac}')
% 
% figure;
% for i=1:X
%     plot(t, y_UKF_MC(3,:,i));hold on;
% end
% title('Courant terminal - UKF')
% xlabel('Temps (sec)')
% ylabel('I_s^2')
% 
% figure;
% for i=1:X
%     plot(t, y_UKF_MC(4,:,i));hold on;
% end
% title('Tension terminale - UKF')
% xlabel('Temps (sec)')
% ylabel('V_s^2')
% 
% figure;
% for i=1:X
%     plot(t, y_UKF_MC(5,:,i));hold on;
% end
% title('Cinquième sortie - UKF')
% xlabel('Temps (sec)')
% ylabel('...')
