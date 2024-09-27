clear

% Code created by Hassan alkhayuon. to compute the stable manifold oa saddel 
% for the Daisy world model. the saddel point computed using COCO
% This work is collaburation with Constantin W. Arnscheidt. 

addpath('../')
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,var)myeventfun(t,var));

%% extracting e_12.
bd1 = coco_bd_read('DaisyW_q1');
ind1_b = coco_bd_idxs(bd1, 'SN');
e12 = coco_bd_col(bd1, 'x');
LL12 = coco_bd_col(bd1, {'L'});

bd3 = coco_bd_read('DaisyW_q3');
e5 = coco_bd_col(bd3, 'x');
LL5 = coco_bd_col(bd3, {'L'});

%% stable manifold of e1

L1_bd = LL12(ind1_b(3):end);
e1_w_bd = e12(1,ind1_b(3):end);
e1_b_bd = e12(2,ind1_b(3):end);

L_start = 1.25;
L_end = 1.53;

res = 250;
L1 = linspace(1.2,1.52,res);
e1_w = interp1(L1_bd,e1_w_bd,L1,"linear");
e1_b = interp1(L1_bd,e1_b_bd,L1,"linear");

figure(10);
hold on
% plot3(L1,e1(1,:),e1(2,:))

del = 1e-5;
tspan = [100 0];
for ind = 1:length(L1)
    odefun = @(t,var)DaisyW(var,L1(ind));
    initcond = [e1_w(ind); e1_b(ind)] + del*[0;1];
    ivpsol = ode45(odefun,tspan,initcond,opts);
    tt = linspace(ivpsol.x(1),ivpsol.x(end),res);

    % plot3(...
    %     L1(ind)*ones(size(ivpsol.y(1,:))), ...
    %     ivpsol.y(1,:),ivpsol.y(2,:))

    ws_temp = deval(ivpsol,tt);
    Ws_e1_w(:,ind) = ws_temp(1,:);
    Ws_e1_b(:,ind) = ws_temp(2,:);
    % plot3(L1(ind).*ones(1,100),ws_temp(1,:),ws_temp(2,:),'k')
end

Ws_e1_surf = surf(L1,Ws_e1_w,Ws_e1_b);
Ws_e1_surf.LineStyle = 'none';
Ws_e1_surf.FaceColor = 'r';
Ws_e1_surf.FaceAlpha = 0.9;

%%

r = 0.1;
par_nonaut = [L_start,L_end,r];
e5_w_start = interp1(LL5,e5(1,:),L_start,"linear");
e5_b_start = interp1(LL5,e5(1,:),L_start,"linear");

initcond = [e5_w_start; e5_b_start];

odefun = @(t,var)DaisyW_nonaut(t,var,par_nonaut);
ivpsol = ode45(odefun,tspan,initcond,opts);

%%
function [check,stop,direction] = myeventfun(t,y)
check = y(1)*y(2)*(y(1) - 0.8)*(y(2) - 0.8);
stop = 1;  % Halt integration
direction = 0;
end