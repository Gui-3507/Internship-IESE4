%% Jacobienne de la fonction de mesure
function H = JacobianMeasFcn(xq, id, iq, vq, rho, Ke, xaf, If, xd)
iex = xaf*If/Ke;
beta = (1-(tanh(rho*Ke))^2)*rho/iex;

H = zeros(5,3);
H(1,1) = -xq*iq;
H(1,2) = vq - xq*id;
H(1,3) = iq;
H(2,1) = -vq;
H(2,2) = -2*xq*iq;
H(2,3) = -id;
H(3,1) = 2*id;
H(3,2) = 2*iq;
H(3,3) = 0;
H(4,1) = 0;
H(4,2) = 2*xq^2*iq;
H(4,3) = 2*vq;
H(5,1) = -beta*xd;
H(5,2) = 0;
H(5,3) = beta;

end

% %% Jacobienne de la fonction de mesure
% function H = JacobianMeasFcn(xq, id, iq, vq, rho, Ke, xaf, If, xd)
% iex = xaf*If/Ke;
% beta = (1-(tanh(rho*Ke))^2)*rho/iex;
% 
% H = zeros(4,3);
% H(1,1) = -xq*iq;
% H(1,2) = vq - xq*id;
% H(1,3) = iq;
% H(2,1) = -vq;
% H(2,2) = -2*xq*iq;
% H(2,3) = -id;
% H(3,1) = 2*id;
% H(3,2) = 2*iq;
% H(3,3) = 0;
% H(4,1) = 0;
% H(4,2) = 2*xq^2*iq;
% H(4,3) = 2*vq;
% % H(5,1) = -beta*xd;
% % H(5,2) = 0;
% % H(5,3) = beta;
% 
% end