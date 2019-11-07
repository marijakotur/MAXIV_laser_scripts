close all
clear all

s=load('30ps_quasi_flattop_spectral.txt');


t=load('30ps_quasi_flattop_temporal.txt');

figure
subplot(3,1,1)
hold on

yyaxis left
plot(s(1,:),s(2,:))
plot(s(1,:),s(3,:),'k--')
ylim([-0.05, 1.05])
ylabel('Intensity [norm.]')

th2 = findall(gca,'Type','line');
for i = 1:length(th2),
   set(th2(i),'linewidth',1.5);
end

yyaxis right
plot(s(1,:),s(4,:)-min(s(4,:)),'-.');
plot(s(1,:),s(5,:)-min(s(5,:)),':');
plot(s(1,:),s(5,:)+s(4,:)-min(s(4,:))-min(s(5,:)),'-');
ylim([-5,45])
xlim([260, 264])
ylabel('Phase [rad]')
th2 = findall(gca,'Type','line');
for i = 1:length(th2),
   set(th2(i),'linewidth',1.5);
end
xlabel('\lambda [nm]')
ylabel('Phase [rad]')
box on
set(gca,'linewidth',1.5);


subplot(3,1,2)
hold on

yyaxis left
plot(t(1,:),t(2,:))
plot(t(1,:),t(3,:),'--')
th2 = findall(gca,'Type','line');
ylim([-0.05, 1.05])
xlim([-inf,inf])
for i = 1:length(th2),
   set(th2(i),'linewidth',1.5);
end
ylabel('Intensity [norm.]')

yyaxis right
plot(t(1,:),t(2,:),'--');
% plot(t(1,:),t(5,:)-min(t(5,:)),':');
plot(t(1,:),t(4,:)-min(t(4,:)),'-');
% ylim([-25,40])
% xlim([259.5, 264.5])
th2 = findall(gca,'Type','line');
for i = 1:length(th2),
   set(th2(i),'linewidth',1.5);
end
xlabel('\lambda [nm]')
ylabel('Phase [rad]')
box on
set(gca,'linewidth',1.5);


saveas(gcf,'pulse_Ew_Et.pdf');