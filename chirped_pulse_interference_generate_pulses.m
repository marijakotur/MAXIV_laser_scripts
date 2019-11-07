function [Et, Et_sh] = generate_pulses(t_fwhm, b, l0, t0, dt, t_span, N)
    %% Generate two gaussian pulse with quadratic spectral phase (i.e. chirp)
%     separted in time by dt.
    
%     t_fwhm: Transform limited pulse duration FWHM in intensity
%     b: quadratic spectral phase parameter (1e-27 range for visible pulses)
%     l0: central wavelength
%     dt: time separation
%     t_span: time span for the generated traces
%     N: number of points, should be enough to resolve the carrier wave
    
    c = 299792458.0;    % Speed of light
    w0 = 2*pi*c/l0;
    tau = t_fwhm/(2*sqrt(log(sqrt(2))));   % Gaussian time constant
    t = linspace(-t_span/2, t_span/2, N);
    t_res = t_span / N;                   % Time resolution
    ph = 0.0
    Eenv = exp(-t.^2/tau^2 + ph) * exp(0j);           % Electric field envelope
    Eenv_sh = exp(-(t+dt).^2/tau^2 + ph) * exp(0j);   % Time shifted e-field envelope
    f_max = 1/(2*t_res);                                   % Maximum frequency for the fft
    w_span = f_max*2*2*pi;
    w = linspace(-w_span/2, w_span/2, N);                 % Generate range of angular frequencies for the fft
    dw = 2*pi*0.441/tau;                                  % Transform limited gaussian pulse w 
    ph = 0.0;
    Eenv_w = fftshift(fft(Eenv));    % FFT to w space
    Eenv_w = Eenv_w.*exp(j*(b*w.^2) + ph);           % Add spectral chirp
    Eenv = ifft(fftshift(Eenv_w));   % iFFT back to t space
    Eenv = Eenv/max(abs(Eenv));
    Et = exp(-1j*w0*t).*Eenv ;                      % Add carrier

    % Same for time shifted pulse
    Eenv_w_sh = fftshift(fft(Eenv_sh));
    Eenv_w_sh = Eenv_w_sh.*exp(1j*(b*w.^2) + ph) ;
    Eenv_sh = ifft(fftshift(Eenv_w_sh));
    Eenv_sh = Eenv_sh/max(abs(Eenv_sh));
    Et_sh = exp(-1j*w0*(t+dt)).*Eenv_sh;
    
end 
