

tbp = 0.44; % for Gaussian pulses
time_duration = 30*1e-15; % in fs
wavelength = 395*1e-9;
c0 = 2.998*1e8;
bandwidth = 0.44*wavelength^2/time_duration/c0;

bw_nm = bandwidth*1e9;
disp([num2str(bw_nm) ' nm'])
