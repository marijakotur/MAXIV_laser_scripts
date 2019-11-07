function Efield = generate_Efield(xdata, x0, sigmax, tdata, t0, sigmat)

% EX = gaussmf(xdata-x0,[2/2.355 0])';
% ET = gaussmf(tdata-t0,[35/2.355 0])';
% Efield = repmat(EX,1,length(tdata)).*repmat(ET,1,length(xdata))';

for ii=1:length(xdata)
    for jj=1:length(tdata)
        Efield(ii,jj) = exp(-(xdata(ii)-x0)^2/(2*sigmax^2)) * exp(-(tdata(jj)-t0)^2/(2*sigmat^2));
    end
end
