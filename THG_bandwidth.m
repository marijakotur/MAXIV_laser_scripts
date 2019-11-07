%bandwidth

lambda = 800; %in nm

% SNLO 800+800=400, type I
delta_k_800 = 17.24 * 10^-7; %in nm^-1 * cm
delta_lambda_800 = lambda^2*delta_k_800/2/pi % in nm per cm
delta_lambda_400 = lambda^2*delta_k_800/2/pi /4;

%thickest crystal
bw_800 = 20;
bw_800/delta_lambda_800

% SNLO 800+400=266, type I