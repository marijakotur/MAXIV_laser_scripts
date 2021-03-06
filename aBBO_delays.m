% a-BBO crystals for pulse stacking
clear all
close all

%% 1mm Kazimir
y1 = load('scandata_2015-11-24_11h20.txt');

%% 2mm Chinese
y2 = load('scandata_2015-11-23_11h50.txt');

% 4mm Chinese
y3 =load('scandata_2015-11-24_13h31.txt');

% 7.5 mm Chinese
y4 =load('scandata_2015-11-24_16h56.txt');

% 15.4 mm Chinese
y5 =load('scandata_2015-11-24_17h14.txt');

figure
hold on
plot((y1(:,1)-.01917)*6666,y1(:,3),'b')
plot((y2(:,1)-.01917)*6666,y2(:,3),'r')
plot((y3(:,1)-.01917)*6666,y3(:,3),'g')
plot((y4(:,1)-.01917)*6666,y4(:,3),'c')
plot((y5(:,1)-.01917)*6666,y5(:,3),'m')



legend('1 mm','2 mm','4 mm','7.5 mm','15.4 mm','Position','best')

figure
plot([15.4,7.5,4,2 ,1],[12.4102    6.0007    3.0292    1.6064    0.6714],'o-')

%% transmission
Theta=[0,  10,  20,  30,  40,  50,  60,  70,  80,  90, 100, 110, 120,...
       130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250,...
       260, 270, 280, 290, 300, 310, 320, 330, 340, 350];
   
   T = [602 551 562; 602 505 510; 602 457 501; 602 384 463; 470 228 303];
   
   figure
   plot([ 1.0000    2.0000    4.0000    7.5000   15.4000],T(:,2)./T(:,1))
 
W16mm=[290, 305, 319, 314, 300, 285, 277, 281, 295, 309, 315, 308, 292,...
       268, 247, 237, 245, 265, 288, 306, 314, 309, 297, 284, 275, 285,...
       295, 314, 317, 310, 294, 272, 250, 242, 248, 267];
 
W4mm=[402, 400, 393, 383, 373, 367, 370, 377, 389, 397, 401, 400, 394,...
       389, 385, 387, 394, 399, 403, 400, 391, 380, 372, 368, 369, 377,...
       390, 398, 402, 400, 393, 388, 385, 386, 391, 401];
 
W0=470;

figure
plot(Theta, W4mm)


fft_length = 512;
W16_padded = [W16mm zeros(fft_length-length(W16mm),1)'];
plot(abs(fft(W16_padded)-mean(W16_padded)))
