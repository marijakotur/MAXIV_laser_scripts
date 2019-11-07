%beam spot size

close all
clear all


% y=imread('10160915-2015-11-05-145850.bmp');
y=imread('10160915-2016-02-11-164254_blue.bmp');
figure
imagesc(y)

figure
hold on
plot(mean(y,2))
plot(mean(y,1))

xdata = 1:length(mean(y,2));
ydata = mean(y,2);
fit=fitcurve_gaussian(xdata',ydata,[68 256 50 0])
% fit=[68 236 50 0];
plot(fit(1) * exp(-(xdata-fit(2)).^2/(2*fit(3)^2)) + fit(4))

%camera sensor 4.8x3.6