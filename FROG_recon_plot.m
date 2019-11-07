%frog reconstruction plot
close all

y=load('frogtrace_2018-10-31_14h51_150mm_12.21_3mm_recon.dat');
t=y(1,:); %in seconds
Et=y(2,:);
phit=y(3,:);
w=y(4,:);
Ew=y(5,:);
phiw=y(6,:);

dw = 2*pi/(max(t)-min(t));
w = linspace(-length(w)/2*dw,length(w)/2*dw,length(w));
lambda = 2*pi*3e8./w 

[t1, Et1, phit1] = array_filter_Intt_phit(t,Et,phit,5000);
[w1, Ew1, phiw1] = array_filter_Intt_phit(w,Ew,phiw,5000);
% w0 = 262;
% w = w+w0;

figure
subplot(2,1,1)
hold on
plotyy(w1,Ew1,w1,phiw1)

subplot(2,1,2)
hold on
plotyy(t1,Et1,t1,phit1)