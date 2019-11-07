tic
%Global Constants ---------------------------------------------------
% All lengths are in mm unless otherwise specified. All times are fs, all freq fs^-1
c = 3.0*10^8;       % speed of light (m/s)
a = 1/(671);        % groove spacing (mm)
lambda = 780e-9;    % central wavelength (m)
f0=c/lambda/1e15;   %central frequency (fs^-1)
theta_i = asin(lambda/(2*(a/1000)));  % incident angle for littrow condition
theta_d = -theta_i;  % diffracted angle for littrow, both in radians
          
p = 1;              % diffraction order
beta = (2*pi*p)/(a*cos(theta_d)*(2*pi*f0)); % in (fs/mm) 1st grating
beta_prime =(2*pi*p)/(a*cos(theta_i)*(2*pi*f0)); % second grating
k0 = (2*pi/780e-6); %central wave vector (rad/mm)
n = 1.5;            %index of the lens
d = 0;              %lens thickness (mm)
fl = 750;           %focal length (mm)
z = fl;             %propagation distance

% Input field in space-time and Fourier domains -------------------------
delta_X = 0.005;              %step size in x (mm)  
delta_T = 1.0;              %step size in t (fs)
fwhm_x = 2.8;
gamma = 1.9;
fwhm_t = 30.5;
tau = (fwhm_t/2)*sqrt(2/log(2));%for field amplitude (fs)
sigma = (fwhm_x/2)*sqrt(2/log(2));%for field amplitude (mm)

N_X = 10.0;                  %spatial extent of x (mm)   
N_T = 50*30*sqrt(2)/2;      %temporal extent of t (fs)

% Field after first grating --------------------------------------------  

[X,T] = meshgrid(-N_X:delta_X:N_X, -N_T:delta_T:N_T);

scale  = 0.0052;
x = [-624:655]*scale;
Ex = load('V:\2006 12 05 A\Edata');
E = interp1(x,Ex,X(1,:));
Ein = kron(exp(-(T(:,1)/tau).^2),E);

clear X
clear T
clear x
clear Ex

%E_lorentz  = sqrt((gamma/2)^2/(X.^2 + (gamma/2)));
% x = X(1,:); %for use as a plotting axis
% E_input = exp(-(x/sigma).^2);% to compare with pulse shaper output
%E_grating1 = exp(-(X/sigma).^2).*exp(-((T-beta*X)/tau).^2);%field after grating
%E_grating1 = sqrt((gamma/2)^2./(X.^2 + (gamma/2))).*exp(-((T-beta*X)/tau).^2);
%E_grating1FT =ifftshift(fft2(E_grating1)); % CHECK ifftshift/fft2!!

clear X
clear T
clear E_grating1

%Propagation from grating 1 to lens 1------------------------------------
fudge=2*pi; % to make beam focus at the focal length - needed to have units of radians
[xi,f]= meshgrid([-1/(2*delta_X):1/2/N_X:1/(2*delta_X)],[-1/(2*delta_T):1/2/N_T:1/(2*delta_T)]);
clear f
E_prop1 = ifft2((E_grating1FT).*exp((i*k0*z)-((i/2/k0)*(z)*(fudge^2*(xi).^2))));
clear xi
clear E_grating1FT


% Action of lens 1 ----------------------------------------------------
[X,T] = meshgrid(-N_X:delta_X:N_X, -N_T:delta_T:N_T);
clear T
E_lens1 =((E_prop1)).*exp(i*k0*n*d - ((i*k0)/2/fl).*(X.^2)); % field just after lens 1
clear E_prop1;
clear X

%Propagation to the mask from lens 1 ------------------------------------
[xi,f]= meshgrid([-1/(2*delta_X):1/2/N_X:1/(2*delta_X)],[-1/(2*delta_T):1/2/N_T:1/(2*delta_T)]);
clear f
E_prop2 = ifft2((fft2(E_lens1)).*exp((i*k0*z)-((i/2/k0)*(z)*(fudge^2*(xi).^2))));%field after propagating to the mask
clear xi;
clear E_lens1;

% Action of the mask --------------------------------------------------

[X,T] = meshgrid(-N_X:delta_X:N_X, -N_T:delta_T:N_T);
clear T
f = 44;
cyclespermm = (1.3*f)/(8.5*4.2); %number of cycles per mm
E_mask = (E_prop2).*exp(i*(pi)*sin(2*pi*cyclespermm*X));%field amplitude after mask
clear E_prop2
clear X


% Propagation after the mask to lens 2 ----------------------------------
[xi,f]= meshgrid([-1/(2*delta_X):1/2/N_X:1/(2*delta_X)],[-1/(2*delta_T):1/2/N_T:1/(2*delta_T)]);
clear f
E_prop3 = ifft2((fft2((E_mask))).*exp((i*k0*z)-((i/2/k0)*(z)*(fudge^2*(xi).^2)))); 
clear E_mask;
clear xi


 % Action of second lens ---------------------------------------------
[X,T] = meshgrid(-N_X:delta_X:N_X, -N_T:delta_T:N_T);
clear T
E_lens2 = ((E_prop3)).*exp(i*k0*n*d - ((i*k0)/2/fl)*(X.^2)); % field just after lens 1
clear E_prop3
clear X

 
% Propagation to second grating from second lens -------------------------
[xi,f]= meshgrid([-1/(2*delta_X):1/2/N_X:1/(2*delta_X)],[-1/(2*delta_T):1/2/N_T:1/(2*delta_T)]);
clear f
E_prop4 = ifft2((fft2((E_lens2))).*exp((i*k0*z)-((i/2/k0)*(z)*(fudge^2*(xi).^2)))); 
clear E_lens2;
clear xi;


% Second Grating
[n,m] = size(E_prop4);
[X,T] = meshgrid(-N_X:delta_X:N_X, -N_T:delta_T:N_T); 
clear T;
for ii=1:m % loop over X
      t_shift = floor(beta_prime*X(1,ii)/delta_T); %How much to shift by
      E_grating2new(:,ii) = circshift(E_prop4(:,ii),[t_shift,0]);
end
clear t_shift
clear E_prop4
clear X


toc
   
% figure
% plot(x,sum((abs(E_grating2new).^2),1)/max(sum((abs(E_grating2new).^2),1)));
% hold on 
% plot(x,((abs(E_input)).^2),'r');
% hold on 
% plot(x,(abs(E_input)).^2/max((abs(E_input)).^2),'r');
% hold on
% plot(x,pulse_before_mask/max(pulse_before_mask),'k')
% xlabel('x(mm)')
% ylabel('Time Integrated Pulse')




