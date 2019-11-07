function Eout = grating_action(Efield, xdata, sigmax, tdata, sigmat, alpha, beta)


% EX = gaussmf(xdata-x0,[2/2.355 0])';
% ET = gaussmf(tdata-t0,[35/2.355 0])';
% Efield = repmat(EX,1,length(tdata)).*repmat(ET,1,length(xdata))';

for ii=1:size(Efield,1)
    Eout(ii,:)=circshift(Efield(ii,:),[0 round(-beta*max(diff(xdata))*xdata(ii))]);
end

%scaling X by alpha
%to do
    
    

% for ii=1:length(xdata)
%     for jj=1:length(tdata)
%         Eout(ii,jj) = exp(-(alpha*xdata(ii))^2/(2*sigmax^2)) * exp(-(tdata(jj)+beta*xdata(ii))^2/(2*sigmat^2));
%     end
% end


end