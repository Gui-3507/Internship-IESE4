%% Fonction de mesure pour les filtres
function y = MeasFcn(iq, id, vq, xq, rho, xaf, Ke, If, xd)
% y(1) : Pac
% y(2) : Qac
% y(3) : Is^2
% y(4) : Vs^2
% y(5) : 

y = zeros(5,1);
iex =  xaf*If/Ke;
y(1) = iq*(vq - xq*id);
y(2) = -xq*iq^2 - id*vq;
y(3) = id^2 + iq^2;
y(4) = xq^2*iq^2 + vq^2;
y(5) = tanh(rho/iex*(vq - xd*id));

end

% %% Fonction de mesure pour les filtres
% function y = MeasFcn(iq, id, vq, xq, rho, xaf, Ke, If, xd)
% % y(1) : Pac
% % y(2) : Qac
% % y(3) : Is^2
% % y(4) : Vs^2
% % y(5) : 
% 
% y = zeros(4,1);
% iex =  xaf*If/Ke;
% y(1) = iq*(vq - xq*id);
% y(2) = -xq*iq^2 - id*vq;
% y(3) = id^2 + iq^2;
% y(4) = xq^2*iq^2 + vq^2;
% 
% 
% end