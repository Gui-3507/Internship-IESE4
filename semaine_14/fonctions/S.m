%% Calcul de la puissance 
function [P, Q] = S(vd, vq, id, iq)

P = vd.*id' + vq.*iq';
Q = vd.*iq' - vq.*id';

end