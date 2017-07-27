    % Fun_leakage_manu.m
% calculates production and consumption leakage of a given tax regime

function [prod_leak, cons_leak, prod_leak_new, cons_leak_new] = Fun_leakage_manu(T)

global sigma beta eta theta gamma jbar Y_rel we we_star pi_c pi_e pi_L pi_cstar pi_Lstar pi_estar G      
tp_prime = T(1); % production tax
tb_prime_frac = T(2); % border tax as a fraction of tp_prime
tb_prime = tb_prime_frac*tp_prime; % border tax adjustment

% I. DEFINE INTERMEDIATE VALUES

% G nails down change in energy price
phat_e = G^((1-beta)/beta);

% change in price of composite energy good
phat_c = phat_e^(1-eta*gamma)*(1+tp_prime)^(1-eta)*(jbar*(1+tp_prime)^(-theta*(1-gamma))...
         + (1-jbar)*(1+tb_prime)^(-theta*(1-gamma)))^(-eta/theta);
phat_cstar = phat_e^(1-eta*gamma)*(jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+(1-jbar))^(-eta/theta);

% price changes enter into changes in consumption shares
pi_cprime = pi_c*(phat_c)^(1-sigma)/(pi_c*(phat_c)^(1-sigma)+1-pi_c);
pi_cstar_prime = pi_cstar*(phat_cstar)^(1-sigma)/(pi_cstar*(phat_cstar)^(1-sigma)+1-pi_cstar);
pihat_c = pi_cprime/pi_c;
pihat_cstar = pi_cstar_prime/pi_cstar;

% change in specialization of manufacturing goods
jbar_prime = jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))/...
             (jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+(1-jbar));

% change in foreign's income (affects home welfare through BTA's)
Yhat_star = pi_Lstar + (1-beta)*pi_estar*phat_e^(1/(1-beta));

% write out change in home income based on parameters and
% quantities that depend on taxes
Yhat = (pi_L + (1-beta)*pi_e*phat_e^(1/(1-beta))+...
       ((tp_prime-tb_prime)/(1+tp_prime))*(1-gamma)*jbar_prime*...
       eta*pi_cstar_prime/Y_rel*Yhat_star)/...
       (1-pi_cprime*(tp_prime/(1+tp_prime)*((1-gamma)*jbar_prime*eta+(1-eta))+...
       tb_prime/(1+tb_prime)*(1-gamma)*(1-jbar_prime)*eta));
   
% II. WRITE OUT CHANGES IN ENERGY USE DIVIDED BY SECTOR AND COUNTRY

Mhat_eh = (((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))/(jbar*...
          ((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+1-jbar))*...
          pihat_c*Yhat/((1+tp_prime)*phat_e);
Mhat_ef = (((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))/(jbar*...
          ((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+1-jbar))*...
          pihat_cstar*Yhat_star/((1+tp_prime)/(1+tb_prime)*phat_e);
Chat_e = pihat_c*Yhat/((1+tp_prime)*phat_e);

Mhat_eh_star = (1/(jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+1-jbar))*...
               pihat_c*Yhat/((1+tb_prime)*phat_e);
Mhat_ef_star = (1/(jbar*((1+tp_prime)/(1+tb_prime))^(-theta*(1-gamma))+1-jbar))*...
               pihat_cstar*Yhat_star/phat_e;
Chat_estar = pihat_cstar*Yhat_star/phat_e;

% III. FORMULATE LEAKAGE

prod_leak = ((1-gamma)*eta*(1-jbar)*(pi_c*(Mhat_eh_star-1)+(pi_cstar/Y_rel*(Mhat_ef_star-1)))+...
            (1-eta)*pi_cstar/Y_rel*(Chat_estar-1))/((1-gamma)*eta*jbar*...
            (pi_c*(1-Mhat_eh)+(pi_cstar/Y_rel*(1-Mhat_ef)))+((1-eta)*pi_c*(1-Chat_e)));
cons_leak = ((1-gamma)*eta*pi_cstar/Y_rel*(jbar*(Mhat_ef-1)+(1-jbar)*(Mhat_ef_star-1))+...
            (1-eta)*pi_cstar/Y_rel*(Chat_estar-1))/((1-gamma)*eta*pi_c*...
            (jbar*(1-Mhat_eh)+(1-jbar)*(1-Mhat_eh_star))+(1-eta)*pi_c*(1-Chat_e));

prod_leak_new = prod_leak/(1-prod_leak);
cons_leak_new = cons_leak/(1-cons_leak);

end

