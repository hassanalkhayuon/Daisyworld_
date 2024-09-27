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


%%  tipping tracking trajectory 
par_nonaut_slow = [L_start; DL; r_slow];
par_nonaut_fast = [L_start; DL; r_fast];
par_nonaut_crit = [L_start; DL; r_crit];

odefun_slow = @(t,var)DaisyW_nonaut(t,var,par_nonaut_slow);
odefun_fast = @(t,var)DaisyW_nonaut(t,var,par_nonaut_fast);
odefun_crit = @(t,var)DaisyW_nonaut(t,var,par_nonaut_crit);

tspan = [-100 100];

basestate = [e5_w(L_start) e5_b(L_start)];

[t_slow,var_slow] = ode45(odefun_slow,tspan,basestate,opts);
[t_fast,var_fast] = ode45(odefun_fast,tspan,basestate,opts);
[t_crit,var_crit] = ode45(odefun_crit,tspan,basestate,opts);

a_w_t_slow = @(tt) interp1(t_slow,var_slow(:,1),tt);
a_b_t_slow = @(tt) interp1(t_slow,var_slow(:,2),tt);

a_w_t_fast = @(tt) interp1(t_fast,var_fast(:,1),tt);
a_b_t_fast = @(tt) interp1(t_fast,var_fast(:,2),tt);

a_w_t_crit = @(tt) interp1(t_crit,var_crit(:,1),tt);
a_b_t_crit = @(tt) interp1(t_crit,var_crit(:,2),tt);


%% slow forcing movie 
end_of_run = length(t_slow);
for ind = 1:end_of_run
    
    L = Lt_slow(t_slow(ind));

    w0 = e0_w(L);
    b0 = e0_b(L);

    w1 = e1_w(L);
    b1 = e1_b(L);

    w2 = e2_w(L);
    b2 = e2_b(L);

    w3 = e3_w(L);
    b3 = e3_b(L);

    w4 = e4_w(L);
    b4 = e4_b(L);

    w5 = e5_w(L);
    b5 = e5_b(L);
    
    time_up_to_now = t_slow(1:ind);

    
    figure(10);
    clf
    set(gcf,'Position',[2202 160 1403 639])

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    cla
    set(gca,'Position',[0.05 0.595 0.44 0.36]);
    set(gca,'FontSize',15);
    hold on
    plot(...
        t_slow,Lt_slow(t_slow),...
        'Color',fainted_green,'LineWidth',2);
    plot(...
        time_up_to_now,Lt_slow(time_up_to_now),...
        'Color',green,'LineWidth',2);
    plot(time_up_to_now(end),Lt_slow(time_up_to_now(end)),...
        '.','MarkerSize',30,'Color',green)
    axis([-100 100 0.7 1.4])
    xticks([-100 -50 0 50 100])
    xticklabels({0 50 100 150})
    yticks([0.8 1 1.2])
    xlabel('$t$','Rotation',0,'Position',[96,0.68,-1])
    ylabel('$L(r t)$','Rotation',0,'Position',[-109,1.3,-1])
    box on

    % second subplot (a_w time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    cla
    set(gca,'Position',[0.05 0.1 0.44 0.36]);
    set(gca,'FontSize',15);
    hold on
    plot(...
        t_slow,e5_w(Lt_slow(t_slow)),...
        'Color',fainted_grey,'LineWidth',2);
    plot(...
        time_up_to_now,a_w_t_slow(time_up_to_now),...
        'Color',green,'LineWidth',2);
    plot(time_up_to_now(end),a_w_t_slow(time_up_to_now(end)),...
        '.','MarkerSize',30,'Color',green)
    axis([-100 100 0 0.8])
    xticks([-100 -50 0 50 100])
    xticklabels({0 50 100 150})
    yticks([0 0.2 0.4 0.6])
    xlabel('$t$','Rotation',0,'Position',[96,-0.015,-1])
    ylabel('$L(r t)$','Rotation',0,'Position',[-109,0.7,-1])
    box on

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    cla
    set(gca,'Position',[0.58 0.1 0.39 0.86]);
    set(gca,'FontSize',15);
    hold on

    if L > 1.194
        plot(w0, b0,'.k','MarkerSize',25)
    else
        plot(w0, b0,'ok','MarkerSize',6,'LineWidth',2)
    end
    % plot(w1, b1,'ok','MarkerSize',6,'LineWidth',2)
    plot(w2, b2,'ok','MarkerSize',6,'LineWidth',2)
    % plot(w3, b3,'ok','MarkerSize',6,'LineWidth',2)
    % plot(w4, b4,'ok','MarkerSize',6,'LineWidth',2)
    plot(w5, b5,'.k','MarkerSize',25)

    plot(...
        a_w_t_slow(time_up_to_now),a_b_t_slow(time_up_to_now),...
        'Color',green,'LineWidth',2);
    plot(a_w_t_slow(time_up_to_now(end)),a_b_t_slow(time_up_to_now(end)),...
        '.','MarkerSize',30,'Color',green)
    


    initcond = [w2+0.001 b2+0.001];
    tspan = [100 0];
    odefun = @(t,var)DaisyW(var,L);

    [t,var_the] = ode45(odefun,tspan,initcond,opts);
    a_w_the = var_the(:,1);
    a_b_the = var_the(:,2);
    axis([0 0.8 0 0.8])
    xticks([0 0.2 0.4 0.6 0.8])
    xticklabels([0 0.2 0.4 0.6])

    yticks([0 0.2 0.4 0.6 0.8])
    yticklabels([0 0.2 0.4 0.6])

    xlabel('$\alpha_w$','Rotation',0,'Position',[0.78,-0.009,-1])
    ylabel('$\alpha_b$','Rotation',0,'Position',[-0.03,0.75,-1])

    box on
    plot(a_w_the,a_b_the,'-k','LineWidth',1.5)
end

%% fast forcing movie 
end_of_run = length(t_fast);
for ind = 1:end_of_run
    
    L = Lt_fast(t_fast(ind));

    w0 = e0_w(L);
    b0 = e0_b(L);

    w1 = e1_w(L);
    b1 = e1_b(L);

    w2 = e2_w(L);
    b2 = e2_b(L);

    w3 = e3_w(L);
    b3 = e3_b(L);

    w4 = e4_w(L);
    b4 = e4_b(L);

    w5 = e5_w(L);
    b5 = e5_b(L);
    
    time_up_to_now = t_fast(1:ind);

    
    figure(10);
    clf
    set(gcf,'Position',[2202 160 1403 639])

    % first subplot (forcing time series)
    h1 = subplot(7,13,[1:6,14:19,27:32]);
    cla
    set(gca,'Position',[0.05 0.595 0.44 0.36]);
    set(gca,'FontSize',15);
    hold on
    plot(...
        t_fast,Lt_fast(t_fast),...
        'Color',fainted_red,'LineWidth',2);
    plot(...
        time_up_to_now,Lt_fast(time_up_to_now),...
        'Color',red,'LineWidth',2);
    plot(time_up_to_now(end),Lt_fast(time_up_to_now(end)),...
        '.','MarkerSize',30,'Color',red)
    axis([-100 100 0.7 1.4])
    xticks([-100 -50 0 50 100])
    xticklabels({0 50 100 150})
    yticks([0.8 1 1.2])
    xlabel('$t$','Rotation',0,'Position',[96,0.68,-1])
    ylabel('$L(r t)$','Rotation',0,'Position',[-109,1.3,-1])
    box on

    % second subplot (a_w time series)
    h2 = subplot(7,13,[53:58,66:71,79:84]);
    cla
    set(gca,'Position',[0.05 0.1 0.44 0.36]);
    set(gca,'FontSize',15);
    hold on
    plot(...
        t_fast,e5_w(Lt_fast(t_fast)),...
        'Color',fainted_grey,'LineWidth',2);
    plot(...
        time_up_to_now,a_w_t_fast(time_up_to_now),...
        'Color',red,'LineWidth',2);
    plot(time_up_to_now(end),a_w_t_fast(time_up_to_now(end)),...
        '.','MarkerSize',30,'Color',red)
    axis([-100 100 0 0.8])
    xticks([-100 -50 0 50 100])
    xticklabels({0 50 100 150})
    yticks([0 0.2 0.4 0.6])
    xlabel('$t$','Rotation',0,'Position',[96,-0.015,-1])
    ylabel('$L(r t)$','Rotation',0,'Position',[-109,0.7,-1])
    box on

    % third subplot (phase space)
    h3 = subplot(7,13,[8:13,21:26,34:39,60:65,73:78,86:91]);
    cla
    set(gca,'Position',[0.58 0.1 0.39 0.86]);
    set(gca,'FontSize',15);
    hold on

    if L > 1.194
        plot(w0, b0,'.k','MarkerSize',25)
    else
        plot(w0, b0,'ok','MarkerSize',6,'LineWidth',2)
    end
    % plot(w1, b1,'ok','MarkerSize',6,'LineWidth',2)
    plot(w2, b2,'ok','MarkerSize',6,'LineWidth',2)
    % plot(w3, b3,'ok','MarkerSize',6,'LineWidth',2)
    % plot(w4, b4,'ok','MarkerSize',6,'LineWidth',2)
    plot(w5, b5,'.k','MarkerSize',25)

    plot(...
        a_w_t_fast(time_up_to_now),a_b_t_fast(time_up_to_now),...
        'Color',red,'LineWidth',2);
    plot(a_w_t_fast(time_up_to_now(end)),a_b_t_fast(time_up_to_now(end)),...
        '.','MarkerSize',30,'Color',red)
    


    initcond = [w2+0.001 b2+0.001];
    tspan = [100 0];
    odefun = @(t,var)DaisyW(var,L);

    [t,var_the] = ode45(odefun,tspan,initcond,opts);
    a_w_the = var_the(:,1);
    a_b_the = var_the(:,2);
    axis([0 0.8 0 0.8])
    xticks([0 0.2 0.4 0.6 0.8])
    xticklabels([0 0.2 0.4 0.6])

    yticks([0 0.2 0.4 0.6 0.8])
    yticklabels([0 0.2 0.4 0.6])

    xlabel('$\alpha_w$','Rotation',0,'Position',[0.78,-0.009,-1])
    ylabel('$\alpha_b$','Rotation',0,'Position',[-0.03,0.75,-1])

    box on
    plot(a_w_the,a_b_the,'-k','LineWidth',1.5)
end
%%
function [check,stop,direction] = myeventfun(t,y)
check = (y(1) - 1)*(y(2) - 1);
stop = 1;  % Halt integration
direction = 0;
end