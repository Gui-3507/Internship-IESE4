%% Fonction d'état pour les filtres
function x_next = StateTransFcn(x)
% x_next(1) : did_dt
% x_next(2) : diq_dt
% x_next(3) : dvq_dt
x_next = zeros(3,1);
x_next(1) = x(1);
x_next(2) = x(2);
x_next(3) = x(3);

end