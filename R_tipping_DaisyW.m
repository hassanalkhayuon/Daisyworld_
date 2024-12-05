clear

% Code created by Hassan alkhayuon. to simulate R-tipping in the 
% nonautonomous the Daisy world model. the base stable point was computed 
% using COCO
% This work is collaburation with Constantin W. Arnscheidt. 

 opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,var)myeventfun(t,var));

 %%

% colour blind friendly colours
red    = [0.86,0.02,0.05];
yellow = [0.87 0.67 0.20];
green  = [0.31,0.70,0.40];
blue   = [0.10,0.40,0.69];
fainted_red   = [0.86,0.02,0.05,0.4];
fainted_green = [0.31,0.70,0.40,0.4];
fainted_blue  = [0.10,0.40,0.69,0.4];
fainted_grey  = [0, 0, 0,0.3];

%% extracting e_12.
bd1 = coco_bd_read('DaisyW_q1');
ind1_b = coco_bd_idxs(bd1, 'SN');
e12 = coco_bd_col(bd1, 'x');
LL12 = coco_bd_col(bd1, {'L'});

bd3 = coco_bd_read('DaisyW_q3');
e5 = coco_bd_col(bd3, 'x');
LL5 = coco_bd_col(bd3, {'L'});


%%

%r = 0.235563107416025;
%r = 0.235563107415577;
r = 1
L_start = 0.8;
DL = 0.5;

par_nonaut = [L_start;DL;r];

e5_w_start = interp1([LL5(15:24), LL5(26:35)],[e5(1,15:24),e5(1,26:35)],L_start);
e5_b_start = interp1([LL5(15:24), LL5(26:35)],[e5(2,15:24),e5(2,26:35)],L_start);

initcond = [0.137556 0.5484426];
% initcond = [0.69036 0];
tspan = [-100 100];

odefun = @(t,var)DaisyW_nonaut(t,var,par_nonaut);
[t,var] = ode45(odefun,tspan,initcond,opts);
a_w = var(:,1);
a_b = var(:,2);
Lt = L_start + (DL./2).*(tanh(r.*t) + 1);

figure(3);
% hold on 
plot(t,a_b,'Color',red)
% xlim([-100 100])
box on 
%%
function [check,stop,direction] = myeventfun(t,y)
check = y(1)*y(2)*(y(1) - 0.8)*(y(2) - 0.8);
stop = 1;  % Halt integration
direction = 0;
end