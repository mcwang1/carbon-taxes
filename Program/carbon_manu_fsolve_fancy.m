    % carbon_manu_fsolve_fancy.m
% Loops over a user-supplied input parameter

% Michael Wang
% August 4, 2016

close all
clear all
clc

% parameters
global sigma beta eta theta gamma jbar Y_rel we we_star pi_c pi_e pi_L pi_cstar pi_Lstar pi_estar G
Opt_para_manu_fancy
Opt_para_manu_derived

% use a cell array to dynamically grab parameters
params = {'sigma';'beta';'eta';'theta';'gamma';'jbar';'Y_rel';'pi_e';'pi_c';'pi_cstar';'G';'tb_prime_frac'};
            
% have user choose the parameter to loop over
[choice,v] = listdlg('PromptString','Select a variable to loop',...
                     'SelectionMode','single','ListString',params);
if(v == 0)
    return;
end

% have user choose values of the parameter to loop over
prompt = {'Starting Value','Step Size','Ending Value'};
param_grid = inputdlg(prompt);

if(isempty(param_grid))
    return;
end

first = str2num(param_grid{1});
step = str2num(param_grid{2});
last = str2num(param_grid{3});

% input grid
param = first:step:last;
n_param = size(param,2);

% results matrices
tp_prime_grid = zeros(1,n_param);
tb_prime_grid = zeros(1,n_param);

% intermediates matrices
phat_e = zeros(1,n_param);
phat_c = zeros(1,n_param);
phat_cstar = zeros(1,n_param);
pi_chat = zeros(1,n_param);
pi_cstar_hat = zeros(1,n_param);
Yhat = zeros(1,n_param);
Yhat_star = zeros(1,n_param);
jbar_hat = zeros(1,n_param);
phat = zeros(1,n_param);
phat_star = zeros(1,n_param);
welfare = zeros(1,n_param);
welfare_star = zeros(1,n_param);

% leakage matrices
prod_leakage = zeros(1,n_param);
cons_leakage = zeros(1,n_param);
prod_leakage_new = zeros(1,n_param);
cons_leakage_new = zeros(1,n_param);

% fsolve parameters
options = optimoptions('fsolve','Display','off','TolFun',1e-6,'MaxFunEvals',1e10,'MaxIter',100);
x0 = .1;

% full or no BTAs?
tb_prime_frac = 0


for i = 1:n_param
    eval([params{choice},'= param(1,i);']);
    Opt_para_manu_derived

    [tp_prime_solution,residuals_goal,flag_goal] = fsolve(@(tp_prime) Fun_goal_manu(tp_prime,tb_prime_frac),x0,options);
    if flag_goal<=0 && sum(sum(abs(residuals_goal)))>1e-6
        fprintf('\n\n Could Not Solve the Problem. Please Try Again.')
        break
    end

    tb_prime_solution = tb_prime_frac*tp_prime_solution;
    tp_prime_grid(1,i) = tp_prime_solution;
    tb_prime_grid(1,i) = tb_prime_solution;

    % intermediates
    [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10] = Opt_intermediates_manu([tp_prime_solution,tb_prime_solution]);
    phat_e(1,i) = i1;
    phat_c(1,i) = i2;
    phat_cstar(1,i) = i3;
    pi_chat(1,i) = i4/pi_c;
    pi_cstar_hat(1,i) = i5/pi_cstar;
    Yhat(1,i) = i6;
    Yhat_star(1,i) = i7;
    jbar_hat(1,i) = i8/jbar;
    phat(1,i) = i9;
    phat_star(1,i) = i10;
    welfare(1,i) = i6/i9;
    welfare_star(1,i) = i7/i10;
    welfare_world(1,i) = wL*welfare(1,i) + wL_star*welfare_star(1,i);

    [l1,l2,l3,l4] = Fun_leakage_manu([tp_prime_solution,tb_prime_frac]);
    prod_leakage(1,i) = l1;
    cons_leakage(1,i) = l2;
    prod_leakage_new(1,i) = l3;
    cons_leakage_new(1,i) = l4;
end

%%

% Figure 1 - results
figure(1)
subplot(2,2,1)
plot(param(1,:),tp_prime_grid(1,:))
title(['Production tax rate vs ' params{choice}])
xlabel(params{choice})
ylabel('tp-prime')

subplot(2,2,2)
plot(param(1,:),prod_leakage(1,:),param(1,:),cons_leakage(1,:))
title(['Leakage vs ' params{choice}])
xlabel(params{choice})
ylabel('leakage')
legend('production leakage','consumption leakage')

subplot(2,2,[3,4])
plot(param(1,:),prod_leakage_new(1,:),param(1,:),cons_leakage_new(1,:))
title(['Modified Leakage vs ' params{choice}])
xlabel(params{choice})
ylabel('modified leakage')
legend('modified prod leak','modified cons leak')

% Figure 2 - Intermediate Values
figure(2)
suptitle(['Intermediate Values vs ' params{choice}])

subplot(2,3,1)
plot(param(1,:),phat_e(1,:),param(1,:),phat(1,:))
hold on
plot(param(1,:),phat_star(1,:))
xlabel(params{choice})
ylabel('price change')
legend('phat_e','phat','phat-star')
hold off

subplot(2,3,2)
plot(param(1,:),phat_c(1,:),param(1,:),phat_cstar(1,:))
xlabel(params{choice})
ylabel('price change')
legend('phat_c','phat_cstar')

subplot(2,3,3)
plot(param(1,:),pi_chat(1,:),param(1,:),pi_cstar_hat(1,:))
xlabel(params{choice})
ylabel('change in consumption share')
legend('pi_chat','pi_chat-star')

subplot(2,3,4)
plot(param(1,:),jbar_hat(1,:))
xlabel(params{choice})
ylabel('jbar-hat')

subplot(2,3,[5 6])
plot(param(1,:),Yhat(1,:),param(1,:),Yhat_star(1,:))
xlabel(params{choice})
ylabel('change in income')
legend('Y-hat','Yhat-star')

figure(3)
plot(param(1,:),welfare(1,:),param(1,:),welfare_world(1,:))
xlabel(params{choice})
ylabel('welfare change')
legend('Home','World')