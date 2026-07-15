%% Filtre de Kalman Étendu
function [x_est_EKF, y_est_EKF] = EKF(Q, R, x, P, rho, Ke, Xd, Xq, Xaf, If, meas, n)

x_est_EKF = zeros(3, n);
y_est_EKF = zeros(5, n);

for i=1:n
    % mise à jour temporelle des états
    x_next = StateTransFcn(x);

    % mise à jour temporelle de la matrice de covariance
    F = JacobianStateFcn();
    P_next = F*P*F' + Q;

    % sortie estimée
    y = MeasFcn(x(2), x(1), x(3), Xq, rho, Xaf, Ke, If, Xd);

    % calcul du gain du filtre
    H = JacobianMeasFcn(Xq, x(1), x(2), x(3), rho, Ke, Xaf, If, Xd);
    K = P_next*H'/(H*P_next*H' + R);

    % mise à jour des états et de la matrice
    x = x_next + K*(meas(:,i) - y);
    P = (eye(3) - K*H)*P_next;
    x_est_EKF(:,i) = x;
    y_est_EKF(:,i) = MeasFcn(x(2), x(1), x(3), Xq, rho, Xaf, Ke, If, Xd);
end

% end

% %% Filtre de Kalman Étendu (avec correction de signe discrète)
% function [x_est_EKF, y_est_EKF] = EKF(Q, R, x, P, rho, Ke, Xd, Xq, Xaf, If, meas, n)
% x_est_EKF = zeros(3, n);
% y_est_EKF = zeros(4, n);
% 
%     for i = 1:n
%         x_next = StateTransFcn(x);
%         F = JacobianStateFcn();
%         P_next = F*P*F' + Q;
% 
%         y = MeasFcn(x(2), x(1), x(3), Xq, rho, Xaf, Ke, If, Xd);
%         H = JacobianMeasFcn(Xq, x(1), x(2), x(3), rho, Ke, Xaf, If, Xd);
%         K = P_next*H'/(H*P_next*H' + R);
% 
%         x_hat = x_next + K*(meas(1:4,i) - y);
%         P_hat = (eye(3) - K*H)*P_next;
% 
%         if (x_hat(3) - Xd*x_hat(1)) < 0   
%             x_hat = -x_hat;                 
%         end
% 
%         x = x_hat;
%         P = P_hat;
% 
%         x_est_EKF(:,i) = x;
%         y_est_EKF(:,i) = MeasFcn(x(2), x(1), x(3), Xq, rho, Xaf, Ke, If, Xd);
%     end
% end

% %% Filtre de Kalman Étendu (avec contrainte inegalité)
% function [x_est_EKF, y_est_EKF] = EKF(Q, R, x, P, rho, Ke, Xd, Xq, Xaf, If, meas, n)
% x_est_EKF = zeros(3, n);
% y_est_EKF = zeros(4, n);   % 4 sorties au lieu de 5 (y5 retirée)
% 
% C = [-Xd, 0, 1];
% N = C';
% 
%     for i = 1:n
%         % mise à jour temporelle des états
%         x_next = StateTransFcn(x);
%         % mise à jour temporelle de la matrice de covariance
%         F = JacobianStateFcn();
%         P_next = F*P*F' + Q;
% 
%         % sortie estimée 
%         y = MeasFcn(x_next(2), x_next(1), x_next(3), Xq, rho, Xaf, Ke, If, Xd);
% 
%         % calcul du gain du filtre 
%         H = JacobianMeasFcn(Xq, x_next(1), x_next(2), x_next(3), rho, Ke, Xaf, If, Xd);
%         K = P_next*H'/(H*P_next*H' + R);
% 
%         % mise à jour des états et de la matrice 
%         x_hat = x_next + K*(meas(1:4,i) - y);
%         P_hat = (eye(3) - K*H)*P_next;
% 
%         if (C*x_hat < 0)
%             % contrainte par projection 
%             D = eye(3) - (N*N')/(N'*N);
%             alpha = 1; % min((-(C*x_next)/(C*D*K*(meas(1:4,i) - y))),1);
%             B = alpha*D;
%             x_proj = x_next + B*K*(meas(1:4,i) - y);
% 
%             G = P_next*H'*K';
%             P_proj = P_next - G*B - B'*G' + B*G*B';
% 
%             % mise à jour 
%             x = x_proj;
%             P = P_proj;
%         else
%             x = x_hat;
% 
%             % P_hat = (P_hat + P_hat') / 2;
%             P = P_hat;
%         end 
% 
%         x_est_EKF(:,i) = x;
%         y_est_EKF(:,i) = MeasFcn(x(2), x(1), x(3), Xq, rho, Xaf, Ke, If, Xd);
%     end
% end