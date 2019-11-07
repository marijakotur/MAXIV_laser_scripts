% Gaussian beam stacking
clear all
close all


t=-10:0.0003:10; %in ps
omega = 7.25*1e3; %in THz

% initial pulse, along dimension 1
t0 = 0;
pulse_fwhm = 2;
sigma_sq = (pulse_fwhm/2.35482).^2;
ENV = exp(-((t-t0).^2/2/sigma_sq));
% ENV2 = exp(-((t-5).^2/2/sigma_sq));
CARRIER = (exp(-i*omega*t));
Et = [ENV.*CARRIER; zeros(1,length(t))]; % electric field [x_comp; y_comp]*carrier*envelope

% crystals
thickness = [4 2 ];
theta = [  45 90 45 90]*pi/180;


for ii=1:length(thickness)
    
    
    % rotation
    ROT =  [cos(theta(ii))  -sin(theta(ii));
        sin(theta(ii))  cos(theta(ii))];
    Et = ROT*Et;
    
    % time shift
    Et = delay_aBBO(t,Et,thickness(ii));
    
    %     ROT =  [cos(theta(ii))  -sin(theta(ii));
    %             sin(theta(ii))  cos(theta(ii))];
    %
    Et = ROT*Et;
    
end

% %plot projection #1
% figure
% hold on
% % plot3(t,zeros(1,length(t)),abs(Et(1,:)).^2,'b--')
% % plot3(t,abs(Et(2,:)).^2,zeros(1,length(t)),'r--')
% % plot3(t,zeros(1,length(t)),abs(Et(1,:)).^2,'b')
% % plot3(t,abs(Et(2,:)).^2,zeros(1,length(t)),'r')
% plot3(t,(Et(1,:).*conj(CARRIER)),(Et(2,:).*conj(CARRIER)),'b')
% view(3)
%
% % change of basis
% theta = 45*pi/180;
% ROT2 = [cos(theta)  -sin(theta);
%         sin(theta)  cos(theta)];
% Et_new = ((ROT*Et)'*ROT^-1)';
% plot3(t,(Et_new(1,:).*conj(CARRIER)),(Et_new(2,:).*conj(CARRIER)),'b')


figure
hold on
fact = sqrt(2)/2;
plot(t,abs(Et(1,:)).^2,'b--')
plot(t,abs(Et(2,:)).^2,'r--')
plot(t,abs(Et(1,:)).^2+abs(Et(2,:)).^2,'k')

show_transmission_data = 0;
if show_transmission_data
    
    %% transmission
    Theta=[0,  10,  20,  30,  40,  50,  60,  70,  80,  90, 100, 110, 120,...
        130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250,...
        260, 270, 280, 290, 300, 310, 320, 330, 340, 350];
    
    W154=[290, 305, 319, 314, 300, 285, 277, 281, 295, 309, 315, 308, 292,...
        268, 247, 237, 245, 265, 288, 306, 314, 309, 297, 284, 275, 285,...
        295, 314, 317, 310, 294, 272, 250, 242, 248, 267]/470;
    
    W75 =[361 358 345 329 311 302 306 316 336 349 360 359 350 340 330 333 ...
        342 354 360 356 346 330 313 303 305 318 337 350 360 360 357 341 330 334 344 353]/470;
    
    W40=[402, 400, 393, 383, 373, 367, 370, 377, 389, 397, 401, 400, 394,...
        389, 385, 387, 394, 399, 403, 400, 391, 380, 372, 368, 369, 377,...
        390, 398, 402, 400, 393, 388, 385, 386, 391, 401]/470;
    
    % W20=[604 601 591 581 571 568 572 580 591 599 602 600 595 590 586 590 593 600 ...
    %     602 597 588 578 566 563 565 570 578 587 590 589 582 574 575 578 587 598]/710;
    
    W20=[341 339 333 332 331 332 333 336 337 339 339 337 336 336 333 338 338 340 ...
        339 337 334 331 329 329 330 330 331 329 327 327 325 327 330 334 337 340]/382;
    
    W10K=[770 771 773 777 777 775 770 765 758 754 756 760 768 772 775 777 775 772 ...
        771 771 773 776 776 774 769 765 759 756 756 759 765 771 775 777 778 774]/851;
    
    W20K = [700 698 708 720 733 741 742 738 730 727 725 732 740 741 741 733 720 707 ...
        698 700 705 717 730 740 740 735 729 726 723 730 736 740 738 729 717 705]/851;
    
    W40K = [619 609 619 641 665 687 695 690 676 662 656 661 674 687 692 684 663 638 ...
        616 609 616 638 664 683 690 687 673 658 653 659 674 686 692 683 663 639]/851;
    
    W80K = [408 396 405 429 460 485 496 491 474 455 443 452 469 489 495 487 464 435 ...
        407 396 405 431 459 485 496 490 474 454 444 453 471 488 496 487 464 433]/692;
    
    figure
    hold on
    f=fitcurve_BRC(Theta,W10K,[0.095 1.5 0.007 1.3 0.8]')
    plot(Theta,W10K,'-o')
    angle=1:360;
    plot(f(1) * sin(2*pi*angle/180+f(2))+ f(3) * sin(4*pi*angle/180+f(4))+f(5))
    
    figure
    subplot(2,1,1)
    hold on
    plot(Theta,W10K,'c-o')
    plot(Theta,W20K,'b-o')
    plot(Theta,W40K,'g-o')
    plot(Theta,W80K,'r-o')
    title('Kasimir')
    legend({'1 mm','2 mm','4 mm','8 mm'},'Location','BestOutside')
    ylim([0.5 1])
    subplot(2,1,2)
    hold on
    plot(Theta,W20,'b-d')
    plot(Theta,W40,'g-d')
    plot(Theta,W75,'r-d')
    plot(Theta,W154,'m-d')
    title('Chinese')
    legend({'2 mm','4 mm','7.5 mm','15.4 mm'},'Location','BestOutside')
end
%
%
%
% %
% % figure
% % hold on
% % polar((Theta-10)*pi/180,W20,'b-o')
% % polar((Theta-10)*pi/180,W20K,'g-o')
% % polar((Theta-10)*pi/180,W40K,'r-o')
% % polar((Theta-10)*pi/180,W80K,'c-o')
% % axis equal