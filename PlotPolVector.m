function []= PlotPolVector(t,Et)

plot3(t,abs(Et(1,:)),zeros(1,length(t)))
hold on
plot3(t,zeros(1,length(t)),abs(Et(1,:)),'r')
fh = gcf;