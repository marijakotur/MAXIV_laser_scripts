%scanning the spectrometer slit to check for spatial chirp
clear all
close all
 
cd('C:\Users\markot\data\2017-11-01\spatial_chirp_scan');
% D = dir('*.txt');
bkgnd_x = 1:100;
bkgnd_y = 1:100;
ind_y = 250:350;
ind_x = 500:900;
d = 0:11;
cm=jet(12);
bkgnd_image = load('background.txt');

for ii = 1:length(d)
%     im = load(D(ii).name);
%     d{ii}=[num2str(str2num(D(ii).name(1:end-2))) ' mm'];
    
%     bkgnd=mean(mean(im(bkgnd_x,bkgnd_y)));
%     im = im - bkgnd;
%     im=im(ind_x,ind_y);
       
    im = load(['hi' num2str(d(ii)) 'mm']);
    im = im - bkgnd_image;
    marg_x = sum(im,2);
    marg_y = sum(im(370:400,475:825),1);
    
          
%     com_x(ii) = sum(marg_x.*(1:length(marg_x))')/sum(marg_x);
    com_y(ii) = sum(marg_y.*(1:length(marg_y)))/sum(marg_y);
%      
%     x_size(ii) = fwhm(1:size(im1,1),marg_x);
    y_size(ii) = fwhm(1:length(marg_y),marg_y);
%      
%     total_signal(ii) = sum(sum(im1));
     
     
%     subplot(2,1,1)
%     hold on
%     plot(marg_x,'color',cm(ii,:))
%     xlabel('pixel')
%     ylabel('x marginal')
% %     legend(d)
    
%     subplot(2,1,2)
    hold on
    plot(marg_y,'color',cm(ii,:))
    axis tight
    xlabel('pixel')
    ylabel('y marginal')
%     legend(d)
     

end

get_data;
ydata = flipud(cell2mat(ydata));
xdata = cell2mat(xdata);


hold on
for ii = 1:length(d)
    fit = fitcurvedemo_gaussian(xdata(ii,:),ydata(ii,:),[max(ydata(ii,:)), 200, 30, 0]);
    com(ii) = fit(2);
    sigma(ii) = fit(3);
    plot(xdata(ii,:),fit(1) * exp(- (xdata(ii,:)-fit(2)).^2/(2*fit(3)^2)) +fit(4),'color',cm(ii,:),'linestyle','--');
end


figure
hold on
plot(com(1:end),'o-')
ylabel('JAI camera pixel')
xlabel('Position in beam [mm]')
title('center of mass')

figure
plot(sigma,'r')