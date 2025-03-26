close all; clear all; clc

%% Import data
V_avoid_full = importdata('V_avoid.mat');
V_clvf_all = importdata('V_clvf.mat');
V_clvf = V_clvf_all(:,:,end);
V_reach = importdata('V_reach.mat');
V_R = flip(V_reach,3);
V_A = flip(V_avoid_full(:,:,1:end),3);
V_reach_avoid = importdata('V_ori_RA.mat');
V_reach2 = importdata('V_reach2.mat');
V_R2 = flip(V_reach2,3);
V_avoid = V_avoid_full(:,:,end);

g_grid = importdata('grid.mat');
obs = V_avoid_full(:,:,1);
goal = V_reach(:,:,1);
goal2 = V_reach2(:,:,1);


% iv_RA_w = importdata('verified_RA_opt_fine.mat');
% iv_RA_c = importdata('verified_RA_opt_fine_correct.mat');

iv_RAS =  importdata('verified_RAS_opt_multi.mat');
% iv_SA = importdata('verified_StabilizeAvoid_fine.mat');

traj = importdata('traj_multi.mat');
traj2 = importdata('traj_multi2.mat');
trajf = importdata('traj_multi_fail1.mat');
trajf2 = importdata('traj_multi_fail2.mat');


%%

figure

subplot(1,3,1)
hold on

% plot(g_grid.xs{1}(ind), g_grid.xs{2}(ind),'r.')
% plot(g_grid.xs{1}(iv_SA), g_grid.xs{2}(iv_SA),'b.')
visSetIm(g_grid,goal,'m',0);
visSetIm(g_grid,V_reach(:,:,end),'m',0);
title('Target and Reach Set')

subplot(1,3,2)
hold on

visSetIm(g_grid,obs,'r',0);
visSetIm(g_grid,V_avoid,'r',0);
title('Obstacle and Avoid Set')


subplot(1,3,3)
hold on
visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0);
% visSetIm(g_grid,V_avoid,'r',0);
visSetIm(g_grid,max(V_avoid,V_reach(:,:,end)),'m',0);
title('RAS and Intersection of Reach and Avoid Set')

%%
fontSize = 25;
titleSize = 20;
exr.LineWidth = 1;

exrval.LineWidth = 1;
exrval.LineStyle = '--' ;

exrClvf.LineWidth = 2;
exrClvf.LineStyle = '--' ;



figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.6]);

hold on
plot(g_grid.xs{1}(iv_RAS), g_grid.xs{2}(iv_RAS),'g.', 'MarkerSize',1)
AVOID = visSetIm(g_grid,V_avoid,'r',0,exrval);
OBS = visSetIm(g_grid,obs,'r',0,exr);
TRAJ = plot(traj(1,:),traj(2,:),'k*');
TRAJ2 = plot(traj2(1,:),traj2(2,:),'k*');
TRAJF = plot(trajf(1,:),trajf(2,:),'r*');
TRAJF2 = plot(trajf2(1,:),trajf2(2,:),'r*');
GOAL = visSetIm(g_grid,goal,'b',0,exr);
GOAL2 = visSetIm(g_grid,goal2,'k',0,exr);

REACH = visSetIm(g_grid,V_reach(:,:,end),'b',0,exrval);
REACH2 = visSetIm(g_grid,V_reach2(:,:,end),'k',0,exrval);

CLVF = visSetIm(g_grid,V_clvf,'m',0.28,exrClvf);
SRCIS = visSetIm(g_grid,V_clvf,'k',-0.01,exrClvf);

% RA = visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0,exrval);
set(gca,'yTick',[-1:1:1]);
set(gca,'xTick',[-1:1:1]);
zx1 = get(gca,'ZTickLabel');
set(gca,'ZTickLabel',zx1,'fontsize',25);
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fontSize );
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fontSize );
title('RAS Set with Two Targets and One Obstacle', 'Interpreter', 'latex', 'FontSize', titleSize );
xlim([-1.5,1.5])
ylim([-1.5,1.5])

lg1 = legend([OBS,GOAL,GOAL2,AVOID,REACH,REACH2,CLVF,SRCIS,TRAJ,TRAJF],...
    {'obstacle', 'target1', 'target2', 'avoid set', 'reach set1',...
    'reach set2','$\mathcal I_M$','$\mathcal I_m$','traj','faild traj'}, ...
    'Interpreter', 'latex', 'FontSize', 18 );


