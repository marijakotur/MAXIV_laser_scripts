clear all
close all

files = dir('uv*.spz');
figure
hold on
for i=1:4
    
    FileName = files(i).name
    fid = fopen(FileName);
    cdata=cell2mat(textscan(fid,'%f %f','delimiter','\t','headerlines',22));
    fclose(fid);
    plot(cdata(:,1),smoothbox(3,cdata(:,2)/max(cdata(:,2))))
    fwhm(cdata(:,1),smoothbox(3,cdata(:,2)/max(cdata(:,2))))
end
leg=legend(files.name);
set(leg,'Interpreter', 'none')
axis([257,267,-inf, inf])
grid on

% files = dir('*3.60*.png');
% figure
% hold on
% wl = load('frogtrace_2018-10-29_16h11_150mm_12.34_2.70mm_wavelengthvector.txt');
% for i=1:length(files)   
%     
%     FileName = files(i).name;
%     data = imread(FileName);
%     data = data - mean(mean(data(end-10:end,end-10,end)));
%     sp = mean(data,1);
%     plot(wl*1e9+2.5,sp/max(sp))
% end


% figure
% hold on
% files = dir('temporal*csv');
% for i=1:length(files)   
%     
%     FileName = files(i).name
%     fid = fopen(FileName);
%     cdata=cell2mat(textscan(fid,'%f %f','delimiter',',','headerlines',1));
%     fclose(fid);
%     plot(cdata(:,1),cdata(:,2)/max(cdata(:,2)))
% end
% legend(files.name)

