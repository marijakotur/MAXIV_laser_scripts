t_span = 10e-12;     % Time span for the generated traces
N = 12000;          % Number of points, should be enough to resolve the carrier wave
t =linspace(-t_span/2, t_span/2, N);


t_fwhm = 80E-15;
b = 17e-27;
l0 = 262e-9;
dt = 0.4e-12;

[Et, Et_sh] = chirped_pulse_interference_generate_pulses(t_fwhm, b, l0,0, dt, t_span, N);
%     t_fwhm: Transform limited pulse duration FWHM in intensity
%     b: quadratic spectral phase parameter (1e-27 range for visible pulses)
%     l0: central wavelength
%     dt: time separation
%     t_span: time span for the generated traces
%     N: number of points, should be enough to resolve the carrier wave


figure(1)
hold on
plot(t*1e12, abs(Et).^2, '-')
plot(t*1e12, abs(Et_sh).^2, '-')
plot(t*1e12, abs(Et+Et_sh).^2)
ax=gca()
% ax.set_xlim(-2e-12, 2e-12)
% ax.grid()
xlabel('Time / ps');
ylabel('Intensity / normalized');
% title('{0:.0f} fs, chirp'.format(t_fwhm*1e15));
grid on

