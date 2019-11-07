%calculate the spot size after a lens

% lambda=input('What was the wavelength in nm? ');
% f=input('What was the focal length in cm? ');
% w_old=input('What was the beam waist in mm? ');
% beam waist is the radius rather then the diameter of the spot

clear all
close all

% all in meters
lambda=262*10^-9;
f=(40:5:50)*10^-2;
w_old=(4:1:10)*10^-3;



for i=1:length(w_old)
    z_old(i)=(pi*w_old(i)^2)/lambda;
    for j=1:length(f)
        factor=sqrt(1+(f(j)^2)./(z_old(i).^2));
        w_new(i,j)=lambda*f(j)/(pi*w_old(i)*factor); %in meters
        w_new_um(i,j)=w_new(i,j)*1e6;
        z(i,j) = pi*w_new(i,j)^2/lambda;
        str{j} = ['f = ' num2str(f(j)*100) ' cm'];
    end
end

figure
subplot(2,1,1)
plot(w_old*1000,w_new_um,'linewidth',2)
legend(str)
xlabel('unfocused spot size [mm]')
ylabel('focused spot size [\mum]')
title(['\lambda=' num2str(lambda*1e9,'%3.0f') 'nm']);
%'f=' num2str(f) 'cm,
axis tight
set(gca,'xminortick','on')
grid on

subplot(2,1,2)
plot(w_new_um,z*1e3,'linewidth',2)
legend(str)
xlabel('focused spot size [um]')
ylabel('Rayleigh range [mm]')

% figure
% surf(f*100,w_old,w_new_um);
% axis tight
% xlabel('focal length [cm]')
% ylabel('spot size before the lens [mm]')
% title('spot size after the lens [\mum]')
% view(0,90)
% colorbar

% figure
% plot(w_old*1000,w_new_um(:,end),'linewidth',3)
%setaxes(gca)



% %calculate intensity
% %P_avg=input('What was the optical power in W? ');
% P_avg=0.001; %1mW
% T=10^-3; %rep rate
% tau=50*10^-15; %pulse duration of 50fs
% P_peak=P_avg*T/tau;
% w_new_cm=w_new*100;
% disp('intensity in W/cm^2 is:')
% %for a gaussian beam the peak intensity on the beam axis
% I=2*P_peak./(pi*w_new_cm.^2);


% plot(w_old,w_new,'linewidth',2);
% xlabel('spot size before lens [mm]')
% ylabel('spor size after lens [um]')

% figure
% plot(w_old,I,'linewidth',2)
% xlabel('w_0 [m]')
% ylabel('peak intensity [W/cm^2')
% title(['int. vs spot size before chamber' num2str(f) ' cm, 50 fs, ' num2str(lambda) ' nm'])
% grid on
%
% figure
% plot(w_new,I,'linewidth',2)
% xlabel('w_0 [m]')
% ylabel('peak intensity [W/cm^2')
% title('int. vs spot size in the chamber' num2str(f) ' cm, 50 fs, ' num2str(lambda) ' nm'])
% grid on


