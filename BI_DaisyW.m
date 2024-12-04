% Code created by Hassan alkhayuon. Tipping diagrams for the 
% Daisyworld model, collaburation with Constantin W. Arnscheidt

clear

%% extracting the base state e5 from COCO bifrucation data.
bd3 = coco_bd_read('DaisyW_q3');
ee5 = coco_bd_col(bd3, 'x');
LL5 = coco_bd_col(bd3, {'L'});

L5 = [LL5(15:24),LL5(26:36)];
e5 = [ee5(:,15:24), ee5(:,26:36)];
%%

% three parameter points 
L_start = 0.8;
DL_max = 1;
r_max = 1;
res = 500;
DL_scan = linspace(0.35,DL_max,res);

r_scan = logspace(-2,1,res);

tol = 1e-3;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10); %,'Events',@(t,var)myeventfun(t,var));

TD_mat = zeros(res,res);

e5_w_start = interp1( L5,e5(1,:),L_start );
e5_b_start = interp1( L5,e5(2,:),L_start );

base_state = [e5_w_start; e5_b_start];
%%
for ind_r = 1:res
    parfor ind_DL = 1:res
        r = r_scan(ind_r);
        DL = DL_scan(ind_DL);
        TD_mat(ind_DL,ind_r) = does_it_tip(L_start,DL,r,base_state);
    end
    disp([ind_r,res])
end
%% computing the threshold 
res_th = 100;
r_th = logspace(-0.9,1,res_th);
DL_init = 0.5;
for ind_th = 1:res_th
    r = r_th(ind_th);
    th_fun = @(DL)does_it_tip(L_start,DL,r,base_state);
    DL_init = fzero(th_fun,DL_init);
    DL_th(ind_th) = DL_init;
    disp(ind_th)
end

%% plotting 
figure(2);
cla
hold on
TD_plot = pcolor(DL_scan, r_scan, TD_mat');
TD_plot.LineStyle = 'none';
TD_plot.FaceAlpha = 0.3;
colormap([1 1 1; 0 0 0])
yscale log
plot(DL_th(1:100),r_th,'Color',[0.6 0.6 0.6])

%% trajectories 
r = 1;
DL = 0.15;

% odefun_mon = @(t,var)DaisyW_monotone(t,var,par_nonaut);
par_nonaut = [L_start; DL; r];

e5_w_start = interp1( L5,e5(1,:),L_start );
e5_b_start = interp1( L5,e5(2,:),L_start );

e5_w_end = interp1( L5,e5(1,:), (L_start + DL) );
e5_b_end = interp1( L5,e5(2,:), (L_start + DL) );

initcond = [e5_w_start; e5_b_start];
T = max( (100/r), 20 );
tspan = [-100 100];

odefun = @(t,var)DaisyW_nonaut(t,var,par_nonaut);
odefun_mon = @(t,var)DaisyW_monotone(t,var,par_nonaut);

[t,var] = ode45(odefun,tspan,initcond,opts);
a_w = var(:,1);
Lt = L_start + (DL./2).*(tanh(r.*t) + 1);

figure(2);
cla
hold on 
plot(t,var(:,2),'k')
figure(3);
cla
hold on 
plot(t,Lt,'k')
%%
function output = does_it_tip(L_start,DL,r,base_state)
tol = 1e-5;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
par_nonaut = [L_start; DL; r];

initcond = base_state;

T = max( (100/r), 1000 );
tspan = [-T, T];
odefun = @(t,var)DaisyW_nonaut(t,var,par_nonaut);

[~,var] = ode45(odefun,tspan,initcond,opts);

output = (var(end,1) <= tol) - (var(end,1) > tol);
end