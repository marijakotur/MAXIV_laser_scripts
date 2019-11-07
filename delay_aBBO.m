function [Et_new] = delay_aBBO(t,Et,thickness)

delay_per_mm = 0.8; %in ps, for aBBO, difference between the axes

%delay
Ew = fft(Et(1,:));
dt=max(diff(t));
w=linspace(-1/dt/2,1/dt/2,length(t));
alpha = 6.29*delay_per_mm/2*thickness; %-6.29 corresponds to a shift of +1ps
Ew_new = abs(Ew).*exp(i*(angle(Ew)+alpha*w));
Et_shifted_1 = ifft(Ew_new);

Ew = fft(Et(2,:));
dt=max(diff(t));
w=linspace(-1/dt/2,1/dt/2,length(t));
alpha = -6.29*delay_per_mm/2*thickness; %-6.29 corresponds to a shift of +1ps
Ew_new = abs(Ew).*exp(i*(angle(Ew)+alpha*w));
Et_shifted_2 = ifft(Ew_new);

Et_new = [Et_shifted_1; Et_shifted_2];