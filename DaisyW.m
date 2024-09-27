function [dvar] = DaisyW(var,par)

% parameters
L = par(1,:);

% daisyworl variables (areas)
a_w = var(1,:);
a_b = var(2,:);

a_g = 1 - a_b - a_w;
    
% albedos
A_b = 0.25;
A_g = 0.5;
A_w = 0.75;
A = A_b.*a_b + A_g.*a_g + A_w.*a_w;

% temperature
sb = 5.67.*10.^(-8); % stefan-boltzmann
S = 917;
q = 0.1.*S./sb;  % heat transfer coefficient (flux/albedo diff) (L neglected)
T_e = (S.*L.*(1-A)./sb).^(0.25);

T_b = (T_e.^4 + q.*(A-A_b)).^0.25;
T_w = (T_e.^4 + q.*(A-A_w)).^0.25;

% rates
b_b = 1 - 0.003265.*(295.5-T_b).^2;
b_w = 1 - 0.003265.*(295.5-T_w).^2;


% b_b = b_b.*(b_b>0);
% b_w = b_w.*(b_w>0);

if b_b<0
    b_b = 0;
end
if b_w<0
    b_w = 0;
end 

g = 0.3;

dvar(1,:) = a_w.*(a_g.*b_w-g);
dvar(2,:) = a_b.*(a_g.*b_b-g);
end