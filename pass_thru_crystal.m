function [new_pulse_parameters] = pass_thru_crystal(time_axis,pulse_parameters,thickness)

thickness 

persistent d
    if isempty(d) %this is only true for the first call to the function
        d=0;
    end
    d = d + 1
%for ii=1:size(old_pulses_parameters,1)

time_shift = 0.8 * thickness(1)/2;

if length(thickness) == 1;
new_pulse_parameters(d) = d;
%  new_pulse_parameters = [pulse_parameters(1)/2 pulse_parameters(2) pulse_parameters(3)+time_shift;
%                          pulse_parameters(1)/2 pulse_parameters(2) pulse_parameters(3)-time_shift];

%  new_pulse_minus = pulse_parameters(1) * ...
%  exp(-((time_axis-time_shift).^2/2/(pulse_parameters(2)/2.35482).^2));
%  new_pulse_plus = pulse_parameters(1) * ...
%  exp(-((time_axis+time_shift).^2/2/(pulse_parameters(2)/2.35482).^2));
%  new_pulses = [new_pulse_minus; new_pulse_plus];
%  size(new_pulses)
else 
  pulse_parameters = [pulse_parameters(1) pulse_parameters(2) pulse_parameters(3)+time_shift]; 
  pass_thru_crystal(time_axis,pulse_parameters,thickness(2:end));
  pulse_parameters = [pulse_parameters(1) pulse_parameters(2) pulse_parameters(3)+time_shift]; 
  pass_thru_crystal(time_axis,pulse_parameters,thickness(2:end));
end