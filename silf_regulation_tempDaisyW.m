clear
% colour blind friendly colours
red    = [0.86,0.02,0.05];
yellow = [0.87 0.67 0.20];
green  = [0.31,0.70,0.40];
blue   = [0.10,0.40,0.69];
fainted_red   = [0.86,0.02,0.05,0.4];
fainted_green = [0.31,0.70,0.40,0.4];
fainted_blue  = [0.10,0.40,0.69,0.4];
fainted_grey  = [0, 0, 0,0.3];

%% bifurcation data 

bd0 = coco_bd_read('DaisyW_q0');
LL0 = coco_bd_col(bd0, {'L'});

% second branch
bd1 = coco_bd_read('DaisyW_q1');
ee1 = coco_bd_col(bd1, 'x');
LL1 = coco_bd_col(bd1, {'L'});

% second branch
bd2 = coco_bd_read('DaisyW_q2');
ee2 = coco_bd_col(bd2, 'x');
LL2 = coco_bd_col(bd2, {'L'});

% third branch
bd3 = coco_bd_read('DaisyW_q3');
ee3 = coco_bd_col(bd3, 'x');
LL3 = coco_bd_col(bd3, {'L'});


L0_1 = LL0(1:8);
L0_2 = LL0(8:20);
L0_3 = LL0(20:end);


L1_1 = LL1(19:32);
L1_2 = LL1(32:50);
L1_3 = LL1(50:68);

a_w1_1 = ee1(1,19:32);
a_w1_2 = ee1(1,32:50);
a_w1_3 = ee1(1,50:68);

a_b1_1 = ee1(2,19:32);
a_b1_2 = ee1(2,32:50);
a_b1_3 = ee1(2,50:68);


L2_1 = LL2(21:33);
L2_2 = LL2(33:53);
L2_3 = LL2(53:70);

a_w2_1 = ee2(1,21:33);
a_w2_2 = ee2(1,33:53);
a_w2_3 = ee2(1,53:70);

a_b2_1 = ee2(2,21:33);
a_b2_2 = ee2(2,33:53);
a_b2_3 = ee2(2,53:70);


L3 = LL3(13:36);
a_w3 = ee3(1,13:36);
a_b3 = ee3(2,13:36);



%%
% import bifurcation data and find L, a_w, a_b

a_g0_1 = ones(size(L0_1));
a_g0_2 = ones(size(L0_2));
a_g0_3 = ones(size(L0_3));

a_g1_1 = 1 - a_b1_1 - a_w1_1;
a_g1_2 = 1 - a_b1_2 - a_w1_2;
a_g1_3 = 1 - a_b1_3 - a_w1_3;

a_g2_1 = 1 - a_b2_1 - a_w2_1;
a_g2_2 = 1 - a_b2_2 - a_w2_2;
a_g2_3 = 1 - a_b2_3 - a_w2_3;

a_g3 = 1 - a_b3 - a_w3;

    
% albedos
A_b = 0.25;
A_g = 0.5;
A_w = 0.75;

A0_1 = A_g.*a_g0_1;
A0_2 = A_g.*a_g0_2;
A0_3 = A_g.*a_g0_3;

A1_1 = A_b.*a_b1_1 + A_g.*a_g1_1 + A_w.*a_w1_1;
A1_2 = A_b.*a_b1_2 + A_g.*a_g1_2 + A_w.*a_w1_2;
A1_3 = A_b.*a_b1_3 + A_g.*a_g1_3 + A_w.*a_w1_3;

A2_1 = A_b.*a_b2_1 + A_g.*a_g2_1 + A_w.*a_w2_1;
A2_2 = A_b.*a_b2_2 + A_g.*a_g2_2 + A_w.*a_w2_2;
A2_3 = A_b.*a_b2_3 + A_g.*a_g2_3 + A_w.*a_w2_3;


A3 = A_b.*a_b3 + A_g.*a_g3 + A_w.*a_w3;



S = 917;

sb = 5.67.*10.^(-8);

T0_1 = (S.*L0_1.*(1-A0_1)./sb).^(1/4);
T0_2 = (S.*L0_2.*(1-A0_2)./sb).^(1/4);
T0_3 = (S.*L0_3.*(1-A0_3)./sb).^(1/4);

T1_1 = (S.*L1_1.*(1-A1_1)./sb).^(1/4);
T1_2 = (S.*L1_2.*(1-A1_2)./sb).^(1/4);
T1_3 = (S.*L1_3.*(1-A1_3)./sb).^(1/4);

T2_1 = (S.*L2_1.*(1-A2_1)./sb).^(1/4);
T2_2 = (S.*L2_2.*(1-A2_2)./sb).^(1/4);
T2_3 = (S.*L2_3.*(1-A2_3)./sb).^(1/4);

T3 = (S.*L3.*(1-A3)./sb).^(1/4);

%% plotting 

figure(2);
cla
hold on

plot(...
    L0_1, T0_1,'-k',...
    L0_2, T0_2,'--k',...
    L0_3, T0_3,'-k','LineWidth',2)

plot(...
    L1_1, T1_1,'--k',...
    L1_2, T1_2,'-k',...
    L1_3, T1_3,'--k','LineWidth',2)

plot(...
    L2_1, T2_1,'--k',...
    L2_2, T2_2,'-k',...
    L2_3, T2_3,'--k','LineWidth',2)

plot(...
    L3, T3,'-k','LineWidth',2)

axis([0.5 1.6 252 337])
box on