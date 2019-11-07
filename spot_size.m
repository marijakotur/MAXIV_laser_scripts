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

% fit=fitcurve_gaussian(length(mean(y,2)),mean(y,2),[100 226 90 8])
% xdata = 1:length(mean(y,2));
% fit=[100 226 90 8];
% plot(fit(1) * exp(- (xdata-fit(2)).^2/(2*fit(3)^2)) + fit(4))