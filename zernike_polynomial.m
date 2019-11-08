%zernike

% rho = 0.1:0.1:1;

x=-1:0.01:1;
y=-1:0.01:1;
n=2;
m=2;

zernike_radial = zeros(length(x));
for i=1:length(x)
    for j=1:length(y)
        rho = sqrt(x(i)^2+y(j)^2);
        phi = atan2(y(i),x(i));
        if rho<1
            if mod(n-m,2) == 0
                for k=0:(n-m)/2
                    zernike_radial(i,j) = zernike_radial(i,j) + (-1)^k*factorial(n-k)/factorial(k)/factorial((n+m)/2-k)/factorial((n-m)/2-k) * rho.^(n-2*k);
                end
            end
        else
            zernike_radial(i,j)=0;
        end
        zernike(i,j) = zernike_radial(i,j) * cos(m*phi);
    end
end

surf(x,y,zernike)
view(2)
axis equal
shading flat
    

