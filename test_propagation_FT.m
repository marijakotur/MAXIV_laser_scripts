clear all
close all

tic
%Global Constants-------------------------------
%All lengths are in mm unless otherwise specified. All times are fs, all freq fs^-1
c=3.0*10^-4;            %speed of light(mm/fs)
a=1/(671);              %groove spacing(mm)
lambda=260e-6;          %central wavelength(mm)
f0=c/lambda;            %central frequency(fs^-1)
thetai=asin(lambda/(2*(a/1000))); %incident angle for littrow condition
thetad=-thetai;         %diffracted angle for littrow, both in radians
p=1;                    %diffractionorder
alpha = 1;              %to correct later for actual deviation from Littrow
beta=(2*pi*p)/(a*cos(thetad)*(2*pi*f0)); %in (fs/mm) 1st grating
betaprime=(2*pi*p)/(a*cos(thetai)*(2*pi*f0)); %second grating
k0=(2*pi/lambda);       %central wavevector(rad/mm)
n=1.5;                  %index of the lens
d=0;                    %lens thickness(mm)
fl=250;                %focal length(mm)
%z=1000;                %propagation distance
fact=1%2*pi;              %needed to have units of radians

%Input field-------------------------
deltaX = 0.04;          %step size in x(mm)
deltaT = 5;             %step size in t(fs)
fwhmx = 2;
gamma = 1.9;
fwhmt = 200;
sigmaT = fwhmt/2.3548;  %for field amplitude(fs)
sigmaX = fwhmx/2.3548;  %for field amplitude(mm)
NX = 9.0;              %spatial extent of x(mm)
NT = 1000;              %temporal extent of t(fs)

X = -NX:deltaX:NX;
T = -NT:deltaT:NT;
T0=0;
X0=0;
xi = -1/(2*deltaX):1/2/NX:1/(2*deltaX);
f = -1/(2*deltaT):1/2/NT:1/(2*deltaT);
k = 2*pi*f/c;

Ein = generate_Efield(X, X0, sigmaX, T, T0, sigmaT);
%Ein = exp(-(T-T0).^2/2*sigmaT^2) .* exp(-(X-X0).^2/2*sigmaX^2);
figure
subplot(2,2,1)
surf(T,X,abs(Ein).^2); shading interp; view(0,90); axis([-500 500 -2 2])
title('Input')

%% Propagation from input to lens1--------------
% EinFT = ifftshift(fft2(Ein));
EinFT = fft2(Ein);
z = 500;
% Eprop1=ifft2((Egrating1FT).*exp((i*k0*z)-((i/2/k0)*(z)*(fact^2*(xi).^2))));
for ii=1:length(xi)
    for jj=1:length(f)
        Eprop1FT(ii,jj) = EinFT(ii,jj)*exp(i*k0*z)*exp(-i*fact^2*xi(ii).^2*z/(2*k0));
    end
end
Eprop1 = ifft2(Eprop1FT);

%% Action of lens 1-------------------------
for ii=1:length(X)
    for jj=1:length(T)
        Elens(ii,jj) = Eprop1(ii,jj)*exp(i*k0*n*d)*exp(-i*k0*X(ii)^2/(2*fl));
    end
end

%% Propagation after lens1--------------
% EinFT = ifftshift(fft2(Ein));
ElensFT = fft2(Elens);

z1 = 250;
Eprop1FT = propagation_action(ElensFT, xi, f, z1, k0);
Eprop1 = ifft2(ElensFT);
subplot(2,2,2)
surf(T,X,abs(Eprop1).^2); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('after lens'); axis([-500 500 -2 2])

z1 = 500;
Eprop2FT = propagation_action(ElensFT, xi, f, z1, k0);
Eprop2 = ifft2(Eprop2FT);
subplot(2,2,3)
surf(T,X,abs(Eprop2).^2); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('after lens'); axis([-500 500 -2 2])

z1 = 1000;
Eprop3FT = propagation_action(ElensFT, xi, f, z1, k0);
Eprop3 = ifft2(Eprop3FT);
subplot(2,2,4)
surf(T,X,abs(Eprop3).^2); shading interp; view(0,90); axis tight; xlabel('Time [fs]'); ylabel('Space [mm]');
title('after lens'); axis([-500 500 -2 2])

toc