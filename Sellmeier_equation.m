function [n, material_str] = Sellmeier_equation(lambda, material_index)

% Sellmeier equation
% lambda must be in micron

if material_index == 1 %a-BBO o
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973&ispreview=yes&tabname=Tutorial
    n = sqrt(2.7471+0.01878./(lambda.^2-0.01822)-0.01354*lambda.^2);
    material_str = 'a-BBO o';
    
elseif material_index == 2 %a-BBO e
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973&ispreview=yes&tabname=Tutorial
    n = sqrt(2.3753+0.01224./(lambda.^2-0.01667)-0.01516*lambda.^2);
    material_str = 'a-BBO e';
    
elseif material_index == 3 % calcite o
    % from http://www.newlightphotonics.com/Birefringent-Crystals/Calcite-Crystals
    n = sqrt(2.69705 + 0.0192064./(lambda.^2-0.01820)-0.0151624*lambda.^2);
    material_str = 'calcite o';
    
elseif material_index == 4 % calcite e
    % from http://www.newlightphotonics.com/Birefringent-Crystals/Calcite-Crystals
    n = sqrt(2.18438 + 0.0087309./(lambda.^2-0.01018)-0.0024411*lambda.^2);
    material_str = 'calcite e';
    
elseif material_index == 5 % fused silica
    %     % from http://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
    %     lambda = lambda*10;
    %     n = sqrt(1+0.6961663*lambda.^2./(lambda.^2-0.06840432) + ...
    %         0.4079426*lambda.^2./(lambda.^2-0.11624142)+0.8974794*lambda.^2./(lambda.^2-9.8961612));
    
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
    n = sqrt(1+0.6962*lambda.^2./(lambda.^2-0.0047) + ...
        0.4079*lambda.^2./(lambda.^2-0.0135)+0.8975*lambda.^2./(lambda.^2-97.934));
    material_str = 'Fused silica';
    
elseif material_index == 6 % Sapphire o
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
    n = sqrt(1+1.4313*lambda.^2./(lambda.^2-0.0047) + ...
        0.4079*lambda.^2./(lambda.^2-0.0143)+5.3414*lambda.^2./(lambda.^2-325.0178));
    material_str = 'Sapphire o';
    
elseif material_index == 7 % Sapphire o
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
    n = sqrt(1+1.5040*lambda.^2./(lambda.^2-0.0055) + ...
        0.5507*lambda.^2./(lambda.^2-0.0148)+6.5927*lambda.^2./(lambda.^2-402.8951));
    material_str = 'Sapphire e';
    
elseif material_index == 8 % CaF2
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
    n = sqrt(1+0.5676*lambda.^2./(lambda.^2-0.0025) + ...
        0.4711*lambda.^2./(lambda.^2-0.0101)+3.8485*lambda.^2./(lambda.^2-1200.5559));
    material_str = 'CaF';
    
elseif material_index == 9 % MgF
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
    n = sqrt(1+0.4134*lambda.^2./(lambda.^2-0.0025) + ...
        0.4711*lambda.^2./(lambda.^2-0.0101)+3.8485*lambda.^2./(lambda.^2-1200.5559));
    material_str = 'MgF';
    
elseif material_index == 10 % MgF
    % from http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
    n = sqrt(1+0.5676*lambda.^2./(lambda.^2-0.0025) + ...
        0.4711*lambda.^2./(lambda.^2-0.0101)+3.8485*lambda.^2./(lambda.^2-1200.5559));
    material_str = 'MgF';
    
elseif material_index == 11 % b-BBO o
    % from http://nathan.instras.com/documentDB/paper-101.pdf
    A = 1.7018379;
    B = 1.0357554;
    C = 0.01800344;
    D = 1.2479989;
    E = 91;
    n = [A + B ./ (1-C./lambda.^2) + D ./ (1 - E./lambda.^2)].^1/2;
    material_str = 'b-BBO o';
    
elseif material_index == 12 % b-BBO e
    % from http://nathan.instras.com/documentDB/paper-101.pdf
    A = 1.5920433;
    B = 0.7816893;
    C = 0.016067891;
    D = 0.8403893;
    E = 91;
    n = [A + B ./ (1-C./lambda.^2) + D ./ (1 - E./lambda.^2)].^1/2;   
    material_str = 'b-BBO e';
    
elseif material_index == 13 % N-BK7
    %http://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=6973
        n = sqrt(1+1.0396*lambda.^2./(lambda.^2-0.006) + ...
        0.2318*lambda.^2./(lambda.^2-0.02)+1.0105*lambda.^2./(lambda.^2-103.5607));
    material_str = 'BK7';
    
else
    print('Material not defined');
end



