% grating angle spread
% pulse shaper

% https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=26
close all
clear all

lambda = (259.4:0.1:264.1)*1e-9;
lambda_0 = 260*1e-9;
N=(3600:50:3850)*10^3;%3846.15
d=1./N;
diff_order = 1;
blaze_wavelength = 240*1E-9;

% figure
% plot(lambda,theta_Littrow)


%gamma = theta_Littrow; %blaze angle

figure(1)
subplot(2,2,1)
hold on
box on

subplot(2,2,2)
hold on

for ii=1:length(d)
    
    % at Littrow theta_in = theta_out = theta_blaze
    theta_blaze(ii) = asin(diff_order*blaze_wavelength./(2*d(ii)));
%     theta_blaze(ii) = blaze_wavelength;
%     theta_in = asin(1*blaze_wavelength./(2*d(ii)));
    theta_in = asin(1*lambda_0./(2*d(ii)))-0.17;
%     theta_in = 25/180*pi;

    
    
    %gamma = -1/2*asin(diff_order*lambda./(2*d(ii)));
    theta_out = asin(-sin(theta_in) + diff_order*lambda/d(ii) ); %for a blazed grating
    
    subplot(2,2,2)
    plot(lambda*1e9,theta_out* 180/pi)
    
    %beam spread at distance L, in mm
    L = 459;
    l(ii) = 2*(tan(max(theta_out)-mean(theta_out)))*L;
    alpha(ii) = (max(theta_out)-min(theta_out));
    
    s{ii} = [num2str(N(ii)*1e-3), 'mm^{-1}'];
    
end

subplot(2,2,1)
plot(N*1e-3,theta_blaze*180/pi,'-o')
size(theta_blaze)
xlabel('N')
ylabel('\theta_{blaze}')

subplot(2,2,2)
box on
xlabel('\lambda [nm]')
ylabel('\theta_{out} [deg]')
axis tight

subplot(2,2,3)
plot(N*1e-3,alpha*180/pi,'bd-')
xlabel('N [mm^{-1}]')
ylabel('Angular spread [deg]')
axis tight
grid on
%legend(s)

subplot(2,2,4)
plot(N*1e-3,l,'-o')
grid on
xlabel('Number of lines')
ylabel('Spatial spread [mm]')


xlabel('N')
