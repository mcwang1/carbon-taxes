% Fun_goal_manu.m
% Computes tp given G, tb

function [residual] = Fun_goal_manu(tp_prime,tb_prime_frac)

global sigma beta eta theta gamma jbar Y_rel we we_star pi_c pi_e pi_L pi_cstar pi_Lstar pi_estar G
tb_prime = tp_prime*tb_prime_frac;

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
       eta*pi_cstar_prime*Yhat_star/Y_rel)/...
       (1-pi_cprime*(tp_prime/(1+tp_prime)*((1-gamma)*jbar_prime*eta+(1-eta))+...
       tb_prime/(1+tb_prime)*(1-gamma)*(1-jbar_prime)*eta));

residual = G - (((1-eta)/(1+tp_prime)+((1-gamma)/(1+tb_prime))*eta*...
      ((1+tb_prime)/(1+tp_prime)*jbar_prime+(1-jbar_prime)))*(we/pi_e)*pi_cprime*Yhat+...
      ((1-eta)+(1-gamma)*eta*((1+tb_prime)/(1+tp_prime)*jbar_prime+(1-jbar_prime)))*...
      (we_star/pi_estar)*pi_cstar_prime*Yhat_star)/phat_e;