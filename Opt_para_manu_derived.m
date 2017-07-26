% Opt_para_manu_derived.m
% Parameters determined by values in Opt_para_manu.m

% back out pi_cstar from market clearing condition
pi_cstar = (pi_e*Y_rel+pi_estar)/(1-eta*gamma) - pi_c*Y_rel
% check for feasibility
if ((pi_cstar < 0) || (pi_cstar > 1))
    disp('Error: unfeasible foreign consumption')
    return
end

% home and foreign's world production shares
we = pi_e*Y_rel/(pi_e*Y_rel + pi_estar)
we_star = pi_estar/(pi_e*Y_rel + pi_estar)

% compute labor shares
pi_L = 1 - (1-beta)*pi_e
pi_Lstar = 1 - (1-beta)*pi_estar
% Check for incomplete specialization of baseline
if ((beta*pi_e + jbar*eta*gamma*pi_c) >= pi_L) || ((beta*pi_estar + (1-jbar)*eta*gamma*pi_cstar) >= pi_Lstar)
    disp('\n\n Outside the bounds of incomplete specialization')
    return
end

% home and foreign's world labor shares
wL = pi_L*Y_rel/(pi_L*Y_rel + pi_Lstar)
wL_star = pi_Lstar/(pi_L*Y_rel + pi_Lstar)