% birefringence
clear all
% close all

% v_g = v_p + k * delta_v_p/delta_k


c=2.998*10^8; %in µm/ps 
step_size = 5*1e-10;  % in m
lambda_nm = (200:step_size:300) * 10^-9; % in µm
lambda = lambda_nm/3;

% Sellmeier equations (lambda in Âµm):

% a-BBO 
% from http://www.newlightphotonics.com/Birefringent-Crystals/alpha-BBO-Crystals
% prob. with a typo
% n_o = sqrt(2.7471+0.01878./(lambda.^2-0.01822)-0.01354*lambda.^2);
% n_e = sqrt(2.37153+0.01224./(lambda.^2-0.01667)-0.01516*lambda.^2);
% mat_str = 'a-BBO'

% from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973&ispreview=yes&tabname=Tutorial
n_o = sqrt(2.7471+0.01878./(lambda.^2-0.01822)-0.01354*lambda.^2);
n_e = sqrt(2.3753+0.01224./(lambda.^2-0.01667)-0.01516*lambda.^2);
mat_str = 'a-BBO';

% %% calcite
% % from http://www.newlightphotonics.com/Birefringent-Crystals/Calcite-Crystals
% n_o = sqrt(2.69705 + 0.0192064./(lambda.^2-0.01820)-0.0151624*lambda.^2);
% n_e = sqrt(2.18438 + 0.0087309./(lambda.^2-0.01018)-0.0024411*lambda.^2);
% mat_str = 'calcite'

subplot(2,1,1)
hold on
plot(lambda,n_o,'b')
plot(lambda,n_e,'r')
legend({[mat_str 'o'],[mat_str 'e']})
xlabel('Wavelength [nm]')
ylabel('Index of refraction')


v_g_o = c./n_o + lambda * c ./n_o.^2 .* diff([n_o(1) n_o])/max(diff(lambda));
v_g_e = c./n_e + lambda * c ./n_e.^2 .* diff([n_e(1) n_e])/max(diff(lambda));

% hold on
% plot(lambda*1000,v_g_o)
% plot(lambda*1000,v_g_e,'r')
% xlabel('\lambda [nm]')

crystal_thickness = 1*1e3; %in µm
delay_o = crystal_thickness./v_g_o;
delay_e = crystal_thickness./v_g_e;

subplot(2,1,2)
hold on
plot(lambda,delay_o*1e3,'b')
plot(lambda,delay_e*1e3,'r')
legend({[mat_str 'GVDo'],[mat_str 'GVDe']})
xlabel('Wavelength [nm]')
ylabel('Delay per mm [fs]')

% figure
% hold on
%plot(lambda,delay_o)
%plot(lambda,delay_e,'r')

figure
plot(lambda,(delay_o-delay_e)*1e3,'r')