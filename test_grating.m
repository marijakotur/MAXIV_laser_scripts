clear all
close all

tic
%Global Constants-------------------------------
%All lengths are in mm unless otherwise specified. All times are fs, all freq fs^-1
c=3.0*10^-4;             %speed of light(mm/fs)
a=1/(600);              %groove spacing(mm)
lambda=800e-9;          %central wavelength(m)
f0=c/lambda/1e3;        %central frequency(fs^-1)
thetai=asin(lambda/(2*(a/1000))); %incident angle for littrow condition
thetad=-thetai;         %diffracted angle for littrow, both in radians
p=1;                    %diffractionorder
alpha = 1;              %to correct later for actual deviation from Littrow
beta=(2*pi*p)/(a*cos(thetad)*2*pi*f0); %in (fs/mm) 1st grating
betaprime=(2*pi*p)/(a*cos(thetai)*2*pi*f0); %second grating
k0=(2*pi/lambda/1000);       %central wavevector(rad/mm)
n=1.5;                  %index of the lens
d=0;                    %lens thickness(mm)
fl=750;                 %focal length(mm)
z=fl;                   %propagation distance
fact=2*pi;              %needed to have units of radians

%Input field-------------------------
deltaX = 0.05;          %step size in x(mm)
deltaT = 4;             %step size in t(fs)
fwhmx = 1;
gamma = 1.9;
fwhmt = 100;
sigmaT = fwhmt/2.3548;  %for field amplitude(fs)
sigmaX = fwhmx/2.3548;  %for field amplitude(mm)
NX = 10.0;              %spatial extent of x(mm)
NT = 1000;              %temporal extent of t(fs)

X = -NX:deltaX:NX;
T = -NT:deltaT:NT;
T0=0;
X0=0;
xi = -1/(2*deltaX):1/2/NX:1/(2*deltaX);
f = -1/(2*deltaT):1/2/NT:1/(2*deltaT);
k = 2*pi*f/c;

%Plotting limits
tmin = -500;
tmax = 500;
smin = -3;
smax = 3;

Ein = generate_Efield(X, X0, sigmaX, T, T0, sigmaT);
Ecurrent = Ein;
clear Ein
%Ein = exp(-(T-T0).^2/2*sigmaT^2) .* exp(-(X-X0).^2/2*sigmaX^2);
figure
subplot(2,4,1)
surf(T,X,(abs(Ecurrent)).^2/max(max(abs(Ecurrent)))); shading interp; view(0,90); axis([tmin tmax smin smax])
title('Input')

%% Field after first grating-------------------
Egrating1 = grating_action(Ecurrent,X, sigmaX,T, sigmaT,alpha, beta);
subplot(2,4,2)
surf(T,X,(abs(Egrating1)).^2/max(max(abs(Egrating1)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 1'); axis([tmin tmax smin smax])
Ecurrent = Egrating1;
clear Egrating1


%% Second Grating
Egrating2 = grating_action(Ecurrent,X, sigmaX,T, sigmaT, 1/alpha, -betaprime/alpha);
subplot(2,4,3)
surf(T,X,(abs(Egrating2)).^2/max(max(abs(Egrating2)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 2'); axis([tmin tmax smin smax])
Ecurrent = Egrating2;


%% Propagation from grating1 to lens1--------------
Egrating1FT = ifftshift(fft2(Ecurrent));
%Eprop1=ifft2((Egrating1FT).*exp((i*k0*z)-((i/2/k0)*(z)*(fact^2*(xi).^2))));
for ii=1:length(xi)
    for jj=1:length(f)
        Eprop1FT(ii,jj) = Egrating1FT(ii,jj)*exp(i*k0*z)*exp(-i*fact^2*xi(ii).^2*z/(2*k0));
    end
end
Eprop1 = ifft2(fftshift(Eprop1FT));
Ecurrent = Eprop1;
clear Eprop1
subplot(2,4,4)
surf(T,X,(abs(Ecurrent)).^2/max(max(abs(Ecurrent)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 1 + prop to lens 1'); axis([tmin tmax smin smax])

%% Propagation from grating1 to lens1--------------
Egrating1FT = ifftshift(fft2(Ecurrent));
%Eprop1=ifft2((Egrating1FT).*exp((i*k0*z)-((i/2/k0)*(z)*(fact^2*(xi).^2))));
z=-z;
for ii=1:length(xi)
    for jj=1:length(f)
        Eprop1FT(ii,jj) = Egrating1FT(ii,jj)*exp(i*k0*z)*exp(-i*fact^2*xi(ii).^2*z/(2*k0));
    end
end
Eprop1 = ifft2(fftshift(Eprop1FT));
Ecurrent = Eprop1;
clear Eprop1
subplot(2,4,5)
surf(T,X,(abs(Ecurrent)).^2/max(max(abs(Ecurrent)))); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('grating 1 + prop to lens 1'); axis([tmin tmax smin smax])
