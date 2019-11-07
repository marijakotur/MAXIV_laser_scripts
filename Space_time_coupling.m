clear all
close all

tic
%Global Constants-------------------------------
%All lengths are in mm unless otherwise specified. All times are fs, all freq fs^-1
c=3.0*10^-4;            %speed of light(mm/fs)
a=1/1500;               %groove spacing(mm)
lambda=260e-6;          %central wavelength(mm)
f0=c/lambda;            %central frequency(fs^-1)
p=1;                    %diffraction order
thetai=asin(p*lambda/(2*a)); %incident angle for Littrow condition
thetad=-thetai;         %diffracted angle for Littrow, both in radians
alpha = 1;              %to correct later for actual deviation from Littrow
beta = p/(a*cos(thetad)*f0); %in (fs/mm) 1st grating
betaprime = p/(a*cos(thetai)*f0); %second grating
k0=2*pi/lambda;       %central wavevector(rad/mm)
n=1.5;                  %index of the lens
d=0;                    %lens thickness(mm)
fl=450;                 %focal length(mm)
z1=450;                  %propagation distance
z2=450;
fact=2*pi;              %needed to have units of radians

%Input field-------------------------
% deltaX = 0.04;        %step size in x(mm)
% deltaT = 10;          %step size in t(fs)
fwhmx = 2;
%gamma = 1.9;
fwhmt = 80;
sigmaT = fwhmt/2.3548;  %for field amplitude(fs)
sigmaX = fwhmx/2.3548;  %for field amplitude(mm)
NX = 2^9;
X_range = 9.0;          %spatial extent of x(mm)
deltaX = 2*X_range / NX;
NT = 2^9;
T_range = 5000;          %temporal extent of t(fs)
deltaT = 2*T_range / NT;

X = linspace(-X_range,X_range,NX);
T = linspace(-T_range,T_range,NT);
T0=0;
X0=0;
xi = linspace(-1/(2*deltaX),1/(2*deltaX),NX);
f = linspace(-1/(2*deltaT),1/(2*deltaT),NT);
k = 2*pi*f/c;



%Plotting limits
tmin = -4000;
tmax = 4000;
smin = -8;
smax = 8;

str = {};

Ein = generate_Efield(X, X0, sigmaX, T, T0, sigmaT);%+generate_Efield(X, X0, sigmaX, T, 300, sigmaT);
Ecurrent = Ein;
%Ein = exp(-(T-T0).^2/2*sigmaT^2) .* exp(-(X-X0).^2/2*sigmaX^2);
figure(1)
subplot(2,4,1)
surf(T,X,(abs(Ecurrent)).^2/max(max(abs(Ecurrent)))); shading interp; view(0,90); axis([tmin tmax smin smax])
title('Input')
figure(2)
hold on
plot(X,sum(abs(Ecurrent).^2,2)/max(sum(abs(Ecurrent).^2,2)),'b')
str{end+1} = 'input'
% 
%% Field after first grating-------------------
Egrating1 = grating_action(Ecurrent,X, sigmaX,T, sigmaT,alpha, beta);
figure(1)
subplot(2,4,2)
surf(T,X,(abs(Egrating1)).^2/max(max(abs(Egrating1)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 1'); axis([tmin tmax smin smax])
Ecurrent = Egrating1;
clear Egrating1

%% Propagation from grating1 to lens1--------------
EcurrentFT = ifftshift(fft2(Ecurrent));
Eprop1FT = propagation_action(EcurrentFT, xi, f, z1, k0);
Eprop1 = ifft2(fftshift(Eprop1FT));
Ecurrent = Eprop1;
clear Eprop1
clear Eprop1FT
figure(2)
plot(X,sum(abs(Ecurrent).^2,2)/max(sum(abs(Ecurrent).^2,2)),'r')
str{end+1} = 'lens1'



%% Action of lens 1-------------------------
for ii=1:length(X)
    for jj=1:length(T)
        Elens(ii,jj) = Ecurrent(ii,jj)*exp(i*k0*n*d)*exp(-i*k0*X(ii)^2/(2*fl));
    end
end
Ecurrent=Elens;
clear Elens
figure(1)
subplot(2,4,3)
surf(T,X,(abs(Ecurrent)).^2/max(max(abs(Ecurrent)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 1 + prop to lens 1'); axis([tmin tmax smin smax])


%% Propagation to the mask from lens 1-------------
EcurrentFT = ifftshift(fft2(Ecurrent));
Eprop2FT = propagation_action(EcurrentFT, xi, f, z1, k0);
Eprop2 = ifft2(fftshift(Eprop2FT));
Ecurrent = Eprop2;
figure(1)
subplot(2,4,4)
surf(T,X,(abs(Eprop2)).^2/max(max(abs(Eprop2)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('lens + prop to mask'); axis([tmin tmax smin smax])
clear Eprop2


%% Action of the mask----------------
%cyclespermm=(1.3*f)/(8.5*4.2);%number of cycles per mm
%Emask = Ecurrent.*exp(i*(pi/2)*sin(2*pi*cyclespermm*X)); % field amplitude after mask

% % blocking a part of the Fourier plane
% for ii=1:size(Ecurrent,1)
%     if and(ii > 0.45*size(Ecurrent,1),ii < 0.55*size(Ecurrent,1))
%         Emask(ii,:)= zeros(size(Ecurrent,2),1);
%     else
%         Emask(ii,:) = Ecurrent(ii,:);
%     end
% end

% for ii=1:size(Ecurrent,1)
%     Emask(ii,:)= Ecurrent(ii,:)*exp(i*0.0005*(ii-size(Ecurrent,1))^2); %chirp
%     %Emask(ii,:)= Ecurrent(ii,:)*exp(-i*cos(ii));
% 
% end

Emask = Ecurrent;
subplot(2,4,5)
surf(T,X,(abs(Emask)).^2/max(max(abs(Emask)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('Mask effect'); axis([tmin tmax smin smax])
figure(2)
plot(X,sum(abs(Ecurrent).^2,2)/max(sum(abs(Ecurrent).^2,2)),'k')
str{end+1} = 'mask'
Ecurrent = Emask;


%% Propagation after the mask to lens2--------------------
EcurrentFT = ifftshift(fft2(Ecurrent));
Eprop3FT = propagation_action(EcurrentFT, xi, f, z2, k0);
Eprop3 = ifft2(Eprop3FT);
Ecurrent = Eprop3;
clear Eprop3


%% Action of second lens-------------------------
for ii=1:length(X)
    for jj=1:length(T)
        Elens2(ii,jj) = Ecurrent(ii,jj)*exp(i*k0*n*d)*exp(-i*k0*X(ii)^2/(2*fl));
    end
end
figure(1)
subplot(2,4,6)
surf(T,X,(abs(Elens2)).^2/max(max(abs(Elens2)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('prop to lens2'); axis([tmin tmax smin smax])
Ecurrent = Elens2;


%% Propagation to second grating from second lens-------
EcurrentFT = (fft2(Ecurrent));
z1 = 450;
Eprop4FT = propagation_action(EcurrentFT, xi, f, z2, k0);
Eprop4 = ifft2(Eprop4FT);
subplot(2,4,7)
surf(T,X,(abs(Eprop4)).^2/max(max(abs(Eprop4)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('prop to grating2'); axis([tmin tmax smin smax])
Ecurrent = Eprop4;
clear Eprop4


%% Second Grating
Egrating2 = grating_action(Ecurrent,X, sigmaX,T, sigmaT, 1/alpha, betaprime/alpha);
subplot(2,4,8)
surf(T,X,(abs(Egrating2)).^2/max(max(abs(Egrating2)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 2'); axis([tmin tmax smin smax])
figure(2)
plot(X,sum(abs(Ecurrent).^2,2)/max(sum(abs(Ecurrent).^2,2)),'g')
str{end+1} = '2nd grating';
legend(str)

toc
% %%
% figure(1)
% subplot(1,2,1)
% hold on
% box on
% plot(T,sum(abs(Ein).^2,1)/max(sum(abs(Ein).^2,1)),'b')
% plot(T,sum(abs(Egrating2).^2,1)/max(sum(abs(Egrating2).^2,1)),'r')
% grid on
% xlabel('Time [fs]'); ylabel('Intensity [arb.]');
% 
% subplot(1,2,2)
% hold on
% box on
% plot(X,sum(abs(Ein).^2,2)/max(sum(abs(Ein).^2,2)),'b')
% plot(X,sum(abs(Egrating2).^2,2)/max(sum(abs(Egrating2).^2,2)),'r')
% grid on
% xlabel('Space [mm]'); ylabel('Intensity [arb.]');

%%
% for ii=1:m %loop over X
%     tshift=floor(betaprime*X(1,ii)/deltaT);%How much to shift by
%     Egrating2(:,ii)=circshift(Eprop4(:,ii),[tshift,0]);
% end

%%
% %% Mode from Pulse Shaper through a lens-----------
% flpinhole=1000;
% %[X,T]=meshgrid(-NX:deltaX:NX,-NT:deltaT:NT);
% %clear T
% Elens=(Egrating2).*exp(i*k0*n*d-((i*k0)/2/flpinhole).*(X.^2));%field just after lens1
% %clear Egrating2;
% %clear X

% %Propagation to pinhole from lens-----------------
% [xi,f]=meshgrid([-1/(2*deltaX):1/2/NX:1/(2*deltaX)],[-1/(2*deltaT):1/2/NT:1/(2*deltaT)]);
% clear f
% Eatpinhole=ifft2((fft2(Elens)).*exp((i*k0*z)-((i/2/k0)*(z)*(fact^2*(xi).^2))));%field after propagating to the mask
% clear xi;
% clear Elens;

% %Field after pinhole------------------------
% [n,m]=size(Eatpinhole);
% Epinhole=[zeros(n,round(m/2)-6) ones(n,11) zeros(n,round(m/2)-6)].*Eatpinhole;
% clear Eatpinhole
% toc

% %% plotting
% figure(1)
% plot(X(:,1),sum((abs(Elens).^2),1)/max(sum((abs(Elens).^2),1)));
% surf(T,X,abs(Elens).^2); shading interp; view(0,90); axis tight
% hold on
% % plot(((abs(Elens)).^2),'r');
% hold on
% plot(X(:,1),(abs(Ein)).^2/max((abs(Ein)).^2),'r');
% hold on
% %plot(X(:,1),pulsebeforemask/max(pulsebeforemask),'k')
% xlabel('x(mm)')
% ylabel('TimeIntegratedPulse')
