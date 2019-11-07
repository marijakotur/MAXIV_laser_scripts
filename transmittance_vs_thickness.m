%bbo transmittance

%thorlabs gives 88% for 6mm

d=1:20;

for ii=1:length(d)
  transmittance(ii) = 0.89.^(d(ii)/6);
end

figure
plot(d,transmittance,'-o')