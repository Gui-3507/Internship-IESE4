%% Filtre de Kalman Unscented
function [x_est_UKF, y_est_UKF] = UKF(x, P, Q, R, Xd, Xq, Xaf, rho, Ke, If, meas, n)
% Définition des constantes et calcul des poids du filtre 
L = length(x);
alpha = 1e-1;
beta = 2;
kappa = 0;
lambda = alpha^2*(L+kappa) - L;

W0_m = lambda/(L + lambda);
W0_c = lambda/(L + lambda) + (1 - alpha^2 + beta);
Wi_m = 1/(2*(L + lambda));
Wi_c = Wi_m;

x_est_UKF = zeros(3, n);
y_est_UKF = zeros(5, n);

for i=1:n
    % calcul des points sigma
    T = chol((L + lambda)*P, "lower");
    chi_ant = [x, x + T, x - T];

    % mise à jour temporelle
    chi_at = zeros(L,2*L + 1);
    for j=1:2*L + 1
        chi_at(:,j) = StateTransFcn(chi_ant(:,j));
    end

    % estimation des états
    % x_pred = chi_at(:,1); %UKO
    x_pred = W0_m*chi_at(:,1);
    for j=2:2*L + 1
        x_pred(:) = x_pred(:) + Wi_m*chi_at(:,j);
    end

    % pré-mise à jour de la matrice de covariance 
    P_pred = Q + W0_c*(chi_at(:,1) - x_pred)*(chi_at(:,1) - x_pred)';
    for j=2:2*L + 1
        P_pred = P_pred + Wi_c*(chi_at(:,j) - x_pred)*(chi_at(:,j) - x_pred)';
    end

    % points sigma dans la fonction de sortie
    gamma = zeros(length(y_est_UKF(:,1)),2*L + 1);
    for j=1:2*L + 1
        gamma(:,j) = MeasFcn(chi_at(2,j), chi_at(1,j), chi_at(3,j), Xq, rho, Xaf, Ke, If, Xd);
    end

    % estimation des sorties
    % y = W0_m*gamma(:,1); % UKO
    y = W0_m*gamma(:,1);
    for j=2:2*L + 1
        y(:) = y(:) + Wi_m*gamma(:,j);
    end

    % mise à jour de la matrice Pyy et Pxy
    Pyy = R + W0_c*(gamma(:,1) - y)*(gamma(:,1) - y)';
    for j=2:2*L + 1
        Pyy = Pyy + Wi_c*(gamma(:,j) - y)*(gamma(:,j) - y)';
    end
    Pxy = W0_c*(chi_at(:,1) - x_pred)*(gamma(:,1) - y)';
    for j=2:2*L + 1
        Pxy = Pxy + Wi_c*(chi_at(:,j) - x_pred)*(gamma(:,j) - y)';
    end

    % gain du filtre
    K = Pxy/Pyy;

    % correction des états estimés et de la matrice de covariance
    x = x_pred + K*(meas(:,i) - y);
    P = P_pred - K*Pyy*K';
    P = nearestSPD(P);

    % sauvegarde des valeurs
    x_est_UKF(:,i) = x;
    y_est_UKF(:,i) = MeasFcn(x(2), x(1), x(3), Xq, rho, Xaf, Ke, If, Xd);
end

end