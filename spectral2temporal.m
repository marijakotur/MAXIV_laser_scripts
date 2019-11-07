function [t,Intt,phit] = spectral2temporal(lambda, absE_lambda, phi_lambda,fftsize)

% lambda -> nu, sample nu uniformly and interpolate
nu = 3e8./lambda*1e9/1e12;
nuI = 3e8/lambda(1)*1e9/1e12;
nuF = 3e8/lambda(end)*1e9/1e12;
d_nu = abs(nuF-nuI)/length(lambda);
nu_uniform = linspace(nuI,nuF,length(lambda));

% omegaIRad=3e8/lambda(1)*1e9/1e12;
% omegaFRad=3e8/lambda(length(lambda))*1e9/1e12;
% dw=(omegaFRad-omegaIRad)/length(lambda);
% omegasRad=omegaIRad:dw:omegaFRad;
% omegaRadPix=3e8./lambda*1e9/1e12;

% absEw=interp1(omegaRadPix, absE_lambda, omegasRad, 'pchip');
% phasew=interp1(omegaRadPix, phi_lambda, omegasRad, 'pchip');
% % absEw, absEl is actually intensity 
% Ew=sqrt(absEw).*exp(i*phasew);

absEnu = interp1(nu, absE_lambda, nu_uniform, 'pchip');
phi_nu=interp1(nu, phi_lambda, nu_uniform, 'pchip');
E_nu=sqrt(absEnu).*exp(i*phi_nu);


pads=zeros(round((fftsize-length(nu))/2),1)';
E_nu_padded=[pads E_nu pads];
size(pads) 
size(E_nu)
% Et_TL=fftshift(fft(fftshift(abs(Ew_padded))));
Et=fftshift(fft(fftshift(E_nu_padded)));
%nu_FFT = d_nu*(1:fftsize);
%d_nu_FFT = 1/fftsize/d_nu;

% dtt_FFT=1/fftsize/d_nu;
t=linspace(-1/2/d_nu,1/2/d_nu,fftsize);
% t=(-1/2/d_nu:dtt_FFT:(1/2/d_nu));

Intt=abs(Et).^2/max(abs(Et).^2);
phit=unwrap(angle(Et));
