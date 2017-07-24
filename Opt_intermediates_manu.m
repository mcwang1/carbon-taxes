% Opt_intermediates_manu.m
% given optimal taxes, calculates intermediate values

function [phat_e,phat_c,phat_cstar,pi_cprime,pi_cstar_prime,...
         Yhat,Yhat_star,jbar_prime,phat,phat_star] = Opt_intermediates_manu(T)

global sigma beta eta theta gamma jbar Y_rel we we_star pi_c pi_e pi_L pi_cstar pi_Lstar pi_estar G      
tp_prime = T(1); % production tax
tb_prime_frac = T(2); % border tax as a fraction of tp_prime
tb_prime = tb_prime_frac*tp_prime; % border tax adjustment
     
% G nails down change in energy price
phat_e = G^((1-beta)/beta);

% change in price of composite energy good
phat_c = phat_e^(1-eta*gamma)*(1+tp_prime)^(1-eta)*(jbar*(1+tp_prime)^(-theta*(1-gamma))...
         + (1-jbar)*(1+tb_prime)^(-theta*(1-gamma)))^(-eta/theta);
phat_cstar = phat_e^(1-eta*gamma)*(jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+(1-jbar))^(-eta/theta);

% price changes enter into changes in consumption shares
pi_cprime = pi_c*(phat_c)^(1-sigma)/(pi_c*(phat_c)^(1-sigma)+1-pi_c);
pi_cstar_prime = pi_cstar*(phat_cstar)^(1-sigma)/(pi_cstar*(phat_cstar)^(1-sigma)+1-pi_cstar);

% change in foreign's income (affects home welfare through BTA's)
Yhat_star = pi_Lstar + (1-beta)*pi_estar*phat_e^(1/(1-beta));

% change in specialization of manufacturing goods
jbar_prime = jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))/...
             (jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+(1-jbar));

% write out change in home income based on parameters and
% quantities that depend on taxes
Yhat = (pi_L + (1-beta)*pi_e*phat_e^(1/(1-beta))+...
       ((tp_prime-tb_prime)/(1+tp_prime))*(1-gamma)*jbar_prime*...
       eta*pi_cstar_prime/Y_rel*Yhat_star)/...
       (1-pi_cprime*(tp_prime/(1+tp_prime)*((1-gamma)*jbar_prime*eta+(1-eta))+...
       tb_prime/(1+tb_prime)*(1-gamma)*(1-jbar_prime)*eta));

% change in home overall price level
phat = (pi_c*(phat_c)^(1-sigma)+1-pi_c)^(1/(1-sigma));

% change in foreign overall price level
phat_star = (pi_cstar*(phat_cstar)^(1-sigma)+1-pi_cstar)^(1/(1-sigma));