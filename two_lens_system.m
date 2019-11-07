%front focal length for a two lens system
clear all
close all

% %+2500mm
% f1 = -152;
% d = 15:0.1:25;
% f2 =160.7;

%-2000mm
f1 = -272.5;
d = 20:0.1:150;
f2 =460.0;


figure(1)
hold on

for i = 1:length(d)
    ffl(i) = f1*(d(i)-f2)./(d(i)-(f1+f2));
 end   
    
    plot(d,ffl,'-o')
    xlabel('f1')
    ylabel('front focal length')

% axis([-inf,inf,0,3500])