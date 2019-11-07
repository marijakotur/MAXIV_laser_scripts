%time bandwidth
%assuming gaussian pulses

tbp = 0.44;

lambda = 392;
lambda_m = lambda * 10^-9;

delta_tau = 25:0.5:100;
delta_tau_s= delta_tau * 10^-15;
c= 3*10^8;

delta_lambda = tbp * lambda_m^2 ./ (c*delta_tau_s) * 10^9;


figure
plot(delta_lambda,delta_tau)
xlabel('\Delta\lambda [nm]')
ylabel('\Delta\tau [fs]')

grid on