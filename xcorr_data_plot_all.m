% plot all cross-correlation scans
clear all
% close all

% data_dir = 'C:\work\xcorr_all\';
% cd(data_dir)
data_files = dir('scandata*14h31.txt'); %_2019-03-18_13

% figure
hold on
ind = 1:length(data_files);
% ind = [1  2  ];
k=0;
for i=1:length(data_files)
    if ismember(i,ind)
        k=k+1;
        y = load(data_files(i).name);
        if min(size(y) ~= 0)
            %         plot(y(:,2)*1e12,y(:,3))
            %         plot(y(:,1)*100,(y(:,3)-mean(y(1:10,3)))/max((y(:,3)-mean(y(1:10,3)))))
%             plot((y(:,1)-y(1,1))*100,k+(y(:,3)-mean(y(1:10,3)))/max((y(:,3)-mean(y(1:10,3)))))
%             plot(y(:,2)*1e12,0*k+(y(:,3)-mean(y(1:10,3)))/max((y(:,3)-mean(y(1:10,3)))))
            z=y(:,3)-mean(y(1:5,3));
            z1=smoothbox(3,z)/max(smoothbox(3,z));
            z1_mean(k) = sum(z1.*y(:,2)*1e12)/sum(z1);
            fit = fitcurve_gaussian(y(:,2)-z1_mean(k),z1,[1 0 3.5 0])
%             subplot(2,1,k)
            plot(y(:,2)*1e12-z1_mean(k),z1+k-1)
            ylabel('Signal [arb]')
            grid on
%             axis([-6, 6, -0.05, 1.05])
            str {k} = data_files(i).name(10:25);
%             legend(data_files(i).name(10:25),'interpreter','none')
        end
    end
end

axis tight
% legend(str)

% xlabel('Stage travel [mm])')
xlabel('Delay [ps])')
% ylabel('Signal [arb]')
box on
%legend('4 mm')
%axis([2 16.5 -inf inf])
% legend(str,'interpreter','none')
grid on
