% Code created by Hassan alkhayuon. Phase portite for the
% Daisyworld model, collaburation with Constantin W. Arnscheidt

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

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,y)myeventfun(t,y));
%% parameter for the nonautonomous system
L_start = 0.8;
DL = 0.5;
r_slow = 0.235563107416025;
         0.235563107416025
r_fast = 0.3;
r_crit = 0.235563107416025;

Lt_slow =@(tt) L_start + (DL./2).*(tanh(r_slow.*tt) + 1);
Lt_fast =@(tt) L_start + (DL./2).*(tanh(r_fast.*tt) + 1);
Lt_crit =@(tt) L_start + (DL./2).*(tanh(r_crit.*tt) + 1);

%% bifurcation data
e0_w =@(LL) 0;
e0_b =@(LL) 0;

% extracting e_1, e_2.
bd1 = coco_bd_read('DaisyW_q1');
ind1_b = coco_bd_idxs(bd1, 'SN');
ee12 = coco_bd_col(bd1, 'x');
LL12 = coco_bd_col(bd1, {'L'});

L1 = LL12(1:37);
e1 = ee12(:,1:37);

e1_w =@(LL) interp1(L1,e1(1,:),LL);
e1_b =@(LL) interp1(L1,e1(2,:),LL);

L2 = LL12(54:end);
e2 = ee12(:,54:end);

e2_w =@(LL) interp1(L2,e2(1,:),LL);
e2_b =@(LL) interp1(L2,e2(2,:),LL);

% extracting e_34.
bd2 = coco_bd_read('DaisyW_q2');
ind2_b = coco_bd_idxs(bd2, 'SN');
ee34 = coco_bd_col(bd2, 'x');
LL34 = coco_bd_col(bd2, {'L'});

L3 = LL34(1:33);
e3 = ee34(:,1:33);

e3_w =@(LL) interp1(L3,e3(1,:),LL);
e3_b =@(LL) interp1(L3,e3(2,:),LL);

L4 = LL34(43:end);
e4 = ee34(:,43:end);

e4_w =@(LL) interp1(L4,e4(1,:),LL);
e4_b =@(LL) interp1(L4,e4(2,:),LL);

% extracting e_5.
bd3 = coco_bd_read('DaisyW_q3');
ind3_b = coco_bd_idxs(bd3, 'BP');
ee5 = coco_bd_col(bd3, 'x');
LL5 = coco_bd_col(bd3, {'L'});

L5 = [LL5(15:24),LL5(26:36)];
e5 = [ee5(:,15:24), ee5(:,26:36)];

e5_w =@(LL) interp1(L5,e5(1,:),LL);
e5_b =@(LL) interp1(L5,e5(2,:),LL);


L = 1.1983;
% L = 0.8;
L = 1.3
figure 
plot(e1_w(L),e1_b(L),'ok','MarkerSize',6,'LineWidth',2)
plot(e2_w(L),e2_b(L),'ok','MarkerSize',6,'LineWidth',2)
plot(e3_w(L),e3_b(L),'ok','MarkerSize',6,'LineWidth',2)
plot(e4_w(L),e4_b(L),'ok','MarkerSize',6,'LineWidth',2)
plot(e5_w(L),e5_b(L),'.k','MarkerSize',25)

%% typical trajectories 
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,y)myeventfun(t,y));
del = 1e-4;
initcond = [e2_w(L); e2_b(L)] + del*[0; 1];
odefun = @(t,var)DaisyW(var,L);
[t,var] = ode45(odefun,[10000 0],initcond,opts);
a_w = var(:,1);
a_b = var(:,2);
figure(3);
hold on 
plot(a_w,a_b,'Color',blue,'LineWidth',2)
area(a_w,a_b,'LineStyle','none')
%%
function [check,stop,direction] = myeventfun(t,y)
check = (y(1) - 10)*(y(2) - 10);
stop = 1;  % Halt integration
direction = 0;
end