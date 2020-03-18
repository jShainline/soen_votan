%% initialize
clc
% clear all
close all

[fontName,fontSize,fontSize_legend,bRGY,scrsz] = f_plotting;
p = f_physicalConstants;

%%
contact_resistance = 30;
m4_resistance = 0.725;
jj1_resistance = 1.45;
num_jjs_vec = [40 80 160];
jc_vec = [2e7 1.5e7 1e7 0.5e7];
d_vec_jja = sqrt(4*40e-6./(pi*jc_vec));
d_vec_jja = linspace(d_vec_jja(1),d_vec_jja(end),8);
d_vec_jjb = sqrt(4*100e-6./(pi*jc_vec));
d_vec_jjb = linspace(d_vec_jjb(1),d_vec_jjb(end),8);
areas_mask_jja = pi*(d_vec_jja/2).^2; 
areas_mask_jjb = pi*(d_vec_jjb/2).^2; 

res_jja_chip29 = {[1302 977 784 651 555 482 429 390],[2490 1881 1503 1241 1053 919 813 733],[4952 3713 2971 2433 2053 1785 1588 1429]};
res_jja_chip28 = {[1352 1018 806 662 562 490 434 390],[2571 1929 1544 1264 1062 917 820 734],[5088 3797 3022 2470 2082 1798 1598 1439]};
res_jja_chip08 = {[1264 948 760 630 538 470 416 377],[2427 1819 1460 1204 1022 889 792 717],[4853 3620 2891 2365 2000 1738 1546 1391]};
res_jja_chip07 = {[1216 910 732 608 519 453 404 367],[2335 1759 1407 1159 986 855 767 694],[4650 3487 2779 2288 1934 1680 1503 1364]};
res_jja_chip01 = {[1132 853 691 575 495 433 388 338],[2179 1648 1321 1103 941 825 742 669],[4360 3287 2634 2183 1856 1619 1450 1310]};

res_jjb_chip30 = {[561 458 382 333 298 273 253 237],[1080 864 725 631 561 511 476 447],[2124 1692 1413 1226 1089 991 915 859]};
res_jjb_chip27 = {[564 456 382 334 300 274 254 239],[1082 869 726 633 562 513 475 447],[2133 1699 1421 1231 1094 996 919 862]};
res_jjb_chip09 = {[536 439 369 325 294 268 251 239],[1034 836 704 616 551 505 471 444],[2044 1633 1373 1197 1068 975 909 858]};
res_jjb_chip06 = {[517 424 361 318 286 263 246 234],[995 809 687 602 542 498 464 440],[1967 1583 1336 1174 1048 961 894 845]};
res_jjb_chip02 = {[485 398 339 299 270 248 232 220],[935 767 650 573 516 474 443 418],[1868 1517 1282 1125 1009 924 862 817]};

for ii = 1:3
    res_jja_chip29{ii}(:) = res_jja_chip29{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jja_chip28{ii}(:) = res_jja_chip28{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jja_chip08{ii}(:) = res_jja_chip08{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jja_chip07{ii}(:) = res_jja_chip07{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jja_chip01{ii}(:) = res_jja_chip01{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    
    res_jjb_chip30{ii}(:) = res_jjb_chip30{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jjb_chip27{ii}(:) = res_jjb_chip27{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jjb_chip09{ii}(:) = res_jjb_chip09{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jjb_chip06{ii}(:) = res_jjb_chip06{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
    res_jjb_chip02{ii}(:) = res_jjb_chip02{ii}(:)-contact_resistance-num_jjs_vec(ii)*(m4_resistance+jj1_resistance);
end

%% converting to Ic and Jc
Delta_Nb = 1.278e-22;
prefactor = pi*Delta_Nb/(2*p.eE);

Ic_jja_chip29 = cell(3,1);
Ic_jja_chip28 = cell(3,1);
Ic_jja_chip08 = cell(3,1);
Ic_jja_chip07 = cell(3,1);
Ic_jja_chip01 = cell(3,1);

Jc_jja_chip29 = cell(3,1);
Jc_jja_chip28 = cell(3,1);
Jc_jja_chip08 = cell(3,1);
Jc_jja_chip07 = cell(3,1);
Jc_jja_chip01 = cell(3,1);

Ic_jjb_chip30 = cell(3,1);
Ic_jjb_chip27 = cell(3,1);
Ic_jjb_chip09 = cell(3,1);
Ic_jjb_chip06 = cell(3,1);
Ic_jjb_chip02 = cell(3,1);

Jc_jjb_chip30 = cell(3,1);
Jc_jjb_chip27 = cell(3,1);
Jc_jjb_chip09 = cell(3,1);
Jc_jjb_chip06 = cell(3,1);
Jc_jjb_chip02 = cell(3,1);

for ii = 1:3
    Ic_jja_chip29{ii}(:) = prefactor./(res_jja_chip29{ii}(:)./num_jjs_vec(ii));
    Ic_jja_chip28{ii}(:) = prefactor./(res_jja_chip28{ii}(:)./num_jjs_vec(ii));
    Ic_jja_chip08{ii}(:) = prefactor./(res_jja_chip08{ii}(:)./num_jjs_vec(ii));
    Ic_jja_chip07{ii}(:) = prefactor./(res_jja_chip07{ii}(:)./num_jjs_vec(ii));
    Ic_jja_chip01{ii}(:) = prefactor./(res_jja_chip01{ii}(:)./num_jjs_vec(ii));
    
    Jc_jja_chip29{ii}(:) = Ic_jja_chip29{ii}(:)./areas_mask_jja(:);
    Jc_jja_chip28{ii}(:) = Ic_jja_chip28{ii}(:)./areas_mask_jja(:);
    Jc_jja_chip08{ii}(:) = Ic_jja_chip08{ii}(:)./areas_mask_jja(:);
    Jc_jja_chip07{ii}(:) = Ic_jja_chip07{ii}(:)./areas_mask_jja(:);
    Jc_jja_chip01{ii}(:) = Ic_jja_chip01{ii}(:)./areas_mask_jja(:);
    
    Ic_jjb_chip30{ii}(:) = prefactor./(res_jjb_chip30{ii}(:)./num_jjs_vec(ii));
    Ic_jjb_chip27{ii}(:) = prefactor./(res_jjb_chip27{ii}(:)./num_jjs_vec(ii));
    Ic_jjb_chip09{ii}(:) = prefactor./(res_jjb_chip09{ii}(:)./num_jjs_vec(ii));
    Ic_jjb_chip06{ii}(:) = prefactor./(res_jjb_chip06{ii}(:)./num_jjs_vec(ii));
    Ic_jjb_chip02{ii}(:) = prefactor./(res_jjb_chip02{ii}(:)./num_jjs_vec(ii));
    
    Jc_jjb_chip30{ii}(:) = Ic_jjb_chip30{ii}(:)./areas_mask_jjb(:);
    Jc_jjb_chip27{ii}(:) = Ic_jjb_chip27{ii}(:)./areas_mask_jjb(:);
    Jc_jjb_chip09{ii}(:) = Ic_jjb_chip09{ii}(:)./areas_mask_jjb(:);
    Jc_jjb_chip06{ii}(:) = Ic_jjb_chip06{ii}(:)./areas_mask_jjb(:);
    Jc_jjb_chip02{ii}(:) = Ic_jjb_chip02{ii}(:)./areas_mask_jjb(:);
end

%% spatial variation

jja_x_coords = [1 1 1 1 1];
jja_y_coords = [1 2 5 6 7];

jjb_x_coords = [2 2 2 2 2];
jjb_y_coords = [1 2 5 6 7];

distance_vec_jja = sqrt( (5*jja_x_coords(:)-0.25).^2 + (5*jja_y_coords(:)-0.25).^2);
distance_vec_jjb = sqrt( (5*jjb_x_coords(:)-0.25).^2 + (5*jjb_y_coords(:)-0.25).^2);

Jc_jja_avg = zeros(5,1);
Jc_jjb_avg = zeros(5,1);
for ii = 1:3
    for jj = 1:8
        Jc_jja_avg(1) = Jc_jja_avg(1)+Jc_jja_chip29{ii}(jj);
        Jc_jja_avg(2) = Jc_jja_avg(2)+Jc_jja_chip28{ii}(jj);
        Jc_jja_avg(3) = Jc_jja_avg(3)+Jc_jja_chip08{ii}(jj);
        Jc_jja_avg(4) = Jc_jja_avg(4)+Jc_jja_chip07{ii}(jj);
        Jc_jja_avg(5) = Jc_jja_avg(5)+Jc_jja_chip01{ii}(jj);
        
        Jc_jjb_avg(1) = Jc_jjb_avg(1)+Jc_jjb_chip30{ii}(jj);
        Jc_jjb_avg(2) = Jc_jjb_avg(2)+Jc_jjb_chip27{ii}(jj);
        Jc_jjb_avg(3) = Jc_jjb_avg(3)+Jc_jjb_chip09{ii}(jj);
        Jc_jjb_avg(4) = Jc_jjb_avg(4)+Jc_jjb_chip06{ii}(jj);
        Jc_jjb_avg(5) = Jc_jjb_avg(5)+Jc_jjb_chip02{ii}(jj);
    end
end
Jc_jja_avg = Jc_jja_avg/24;
Jc_jjb_avg = Jc_jjb_avg/24;

%% total resistance plots

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jja*1e12,res_jja_chip29{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jja*1e12,res_jja_chip29{2},'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jja*1e12,res_jja_chip29{3},'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('Room Temp. Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jja, Chip 29','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb*1e12,res_jjb_chip30{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb*1e12,res_jjb_chip30{2},'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb*1e12,res_jjb_chip30{3},'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('Room Temp. Resistance [\Omega]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jjb, Chip 30','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

% k1 = gtext(figureCaptions(1:length(figureCaptions)));
% set(k1,'FontSize',fontSize_legend,'FontName','Times')

%% resistance per jj plots

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jja*1e12,res_jja_chip29{1}/40,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jja*1e12,res_jja_chip29{2}/80,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jja*1e12,res_jja_chip29{3}/160,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('Room Temp. Resistance per JJ [\Omega]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jja, Chip 29','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb*1e12,res_jjb_chip30{1}/40,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb*1e12,res_jjb_chip30{2}/80,'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb*1e12,res_jjb_chip30{3}/160,'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('Room Temp. Resistance per JJ [\Omega]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jjb, Chip 30','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

%% Ic plots

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jja*1e12,1e6*Ic_jja_chip29{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jja*1e12,1e6*Ic_jja_chip29{2},'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jja*1e12,1e6*Ic_jja_chip29{3},'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('I_c [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jja, Chip 29','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb*1e12,1e6*Ic_jjb_chip30{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb*1e12,1e6*Ic_jjb_chip30{2},'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb*1e12,1e6*Ic_jjb_chip30{3},'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('I_c [\mu A]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jjb, Chip 30','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

%% Jc plots

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jja*1e12,1e-7*Jc_jja_chip29{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jja*1e12,1e-7*Jc_jja_chip29{2},'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jja*1e12,1e-7*Jc_jja_chip29{3},'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('J_c [kA/cm^2]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jja, Chip 29','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(areas_mask_jjb*1e12,1e-7*Jc_jjb_chip30{1},'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
hold on
plot(areas_mask_jjb*1e12,1e-7*Jc_jjb_chip30{2},'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
plot(areas_mask_jjb*1e12,1e-7*Jc_jjb_chip30{3},'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
legend('40 jjs','80 jjs','160 jjs')
ylabel('J_c [kA/cm^2]','FontSize',fontSize,'FontName','Times')
xlabel('Junction Area [\mu m^2]','FontSize',fontSize,'FontName','Times')
title('jjb, Chip 30','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

%% spatial variation

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(distance_vec_jja,1e-7*Jc_jja_avg,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
ylabel('Average jj_a J_c [kA/cm^2]','FontSize',fontSize,'FontName','Times')
xlabel('Radial Distance [mm]','FontSize',fontSize,'FontName','Times')
% title('Chip 29','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
plot(distance_vec_jjb,1e-7*Jc_jjb_avg,'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
ylabel('Average jj_b J_c [kA/cm^2]','FontSize',fontSize,'FontName','Times')
xlabel('Radial Distance [mm]','FontSize',fontSize,'FontName','Times')
% title('Chip 29','FontSize',fontSize,'FontName','Times')
set(gca,'FontSize',fontSize,'FontName',fontName)

% % k1 = gtext(figureCaptions(1:length(figureCaptions)));
% % set(k1,'FontSize',fontSize_legend,'FontName','Times')