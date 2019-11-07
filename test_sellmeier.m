%test Sellmeier function
close all
clear all

lambda = 255:1:850;
figure
hold on
num =13;
cm = hsv(num);
k=1
ind = [5,8]
for ii=1:length(ind)
    [n mat_str]=Sellmeier_equation(lambda/1000,ind(ii));
    plot(lambda,n,'color',cm(ii,:))
    str{k} = mat_str;
    k=k+1;   
end
    
legend(str)
grid on

