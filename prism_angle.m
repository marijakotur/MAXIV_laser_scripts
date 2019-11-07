% fused silica prism
close all
clear all

lambda = (255:0.1:265)/1000; %in micron

% %fused silica
n = sqrt(1 + 0.6961663*lambda.^2./(lambda.^2-0.0684043^2) + 0.4079426*lambda.^2./(lambda.^2-0.1162414^2) + 0.8974794*lambda.^2./(lambda.^2-9.896161^2));

%sapphire
% n = sqrt(1 + 1.4313493*lambda.^2./(lambda.^2-0.0726631^2) + 0.65054713*lambda.^2./(lambda.^2-0.1193242^2) + 5.3414021*lambda.^2./(lambda.^2-18.028251^2));

n_air = 1.000293;


%Brewster angle
theta_B = atan(n/n_air);


%angle of minimum deviation

% a ray coming in at theta_B refracts at:
theta_B_prime = asin(n_air./n .* sin(theta_B));
A = 2*theta_B_prime;

% minimum deviation
% A=60*pi/180; %apex angle
d = 2 * (-A/2 + asin(n.*sin(A/2))); %minimum deviation angle


hold on
plot(lambda,theta_B*180/pi)
plot(lambda,A*180/pi,'r')
legend('Brewster','Apex')
xlabel('Wavelength')
ylabel('Angle')
grid on

% plot(lambda,d*180/pi)


