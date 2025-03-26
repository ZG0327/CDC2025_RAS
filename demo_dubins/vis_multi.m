close all; clear all; clc

%% Import data
V_avoid_full = importdata('V_avoid.mat');
V_clvf_all = importdata('V_clvf.mat');
V_clvf = V_clvf_all(:,:,:,end);
V_reach = importdata('V_reach.mat');
V_R = flip(V_reach,3);
V_A = flip(V_avoid_full(:,:,:,1:end),3);
V_reach_avoid = importdata('V_ori_RA.mat');
V_reach2 = importdata('V_reach2.mat');
V_R2 = flip(V_reach2,3);
V_avoid = V_avoid_full(:,:,:,end);

g_grid = importdata('grid.mat');
obs = V_avoid_full(:,:,:,1);
goal = V_reach(:,:,:,1);
goal2 = V_reach2(:,:,:,1);

[ g_p, goal_p] = proj(g_grid,goal,[0,0,1],'min');
[~, goal2_p] = proj(g_grid,goal2,[0,0,1],'min');
[~, obs_p] = proj(g_grid,obs,[0,0,1],'min');
[~, clvf_p] = proj(g_grid,V_clvf,[0,0,1],'min');
[~, ra] = proj(g_grid,V_reach_avoid,[0,0,1],'min');

% iv_RA_w = importdata('verified_RA_opt_fine.mat');
% iv_RA_c = importdata('verified_RA_opt_fine_correct.mat');

iv_RAS =  importdata('verified_RAS_opt_multi.mat');
% iv_SA = importdata('verified_StabilizeAvoid_fine.mat');

traj = importdata('traj_multi.mat');
traj2 = importdata('traj_multi2.mat');

%%

% figure
% 
% subplot(1,3,1)
% hold on
% 
% % plot(g_grid.xs{1}(ind), g_grid.xs{2}(ind),'r.')
% % plot(g_grid.xs{1}(iv_SA), g_grid.xs{2}(iv_SA),'b.')
% visSetIm(g_grid,goal,'m',0);
% visSetIm(g_grid,V_reach(:,:,end),'m',0);
% title('Target and Reach Set')
% 
% subplot(1,3,2)
% hold on
% 
% visSetIm(g_grid,obs,'r',0);
% visSetIm(g_grid,V_avoid,'r',0);
% title('Obstacle and Avoid Set')
% 
% 
% subplot(1,3,3)
% hold on
% visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0);
% % visSetIm(g_grid,V_avoid,'r',0);
% visSetIm(g_grid,max(V_avoid,V_reach(:,:,end)),'m',0);
% title('RAS and Intersection of Reach and Avoid Set')

%%
fontSize = 25;
titleSize = 20;
exr.LineWidth = 1;

exrval.LineWidth = 1;
exrval.LineStyle = '-.' ;
exrval.applyLight = 1;

exrClvf.LineWidth = 2;
exrClvf.LineStyle = '--' ;

figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.6]);

subplot(1,2,1)
set(gca,'unit','normalized','position',[0.1,0.3,0.35,0.6])

hold on
plot3(g_grid.xs{1}(iv_RAS), g_grid.xs{2}(iv_RAS),g_grid.xs{3}(iv_RAS),'g.', 'MarkerSize',3)
% AVOID = visSetIm(g_grid,V_avoid,'r',0,exrval);
OBS = visSetIm(g_grid,obs,'r',0,exr);
TRAJ = plot3(traj(1,:),traj(2,:),traj(3,:),'k*')
plot3(traj2(1,1:152),traj2(2,1:152),traj2(3,1:152),'k*')
GOAL = visSetIm(g_grid,goal,'b',0,exr);
GOAL2 = visSetIm(g_grid,goal2,'k',0,exr);

% REACH = visSetIm(g_grid,V_reach(:,:,:,end),'b',0,exrval);
% REACH2 = visSetIm(g_grid,V_reach2(:,:,:,end),'k',0,exrval);

CLVF = visSetIm(g_grid,V_clvf,'m',0.1,exrClvf);

% RA = visSetIm(g_grid,V_reach_avoid(:,:,end),'k',0,exrval);
set(gca,'yTick',[-3:4:5]);
set(gca,'xTick',[-5:5:5]);
set(gca,'zTick',[-3:3:3]);

zx1 = get(gca,'ZTickLabel');
set(gca,'ZTickLabel',zx1,'fontsize',25);
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', fontSize );
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', fontSize );
zlabel('$x_3$', 'Interpreter', 'latex', 'FontSize', fontSize );

sgtitle('RAS Set with Two Targets and One Obstacle', 'Interpreter', 'latex', 'FontSize', titleSize );
xlim([-5,5]);
ylim([-3,5]);
zlim([-pi,pi]);



subplot(1,2,2)
set(gca,'unit','normalized','position',[0.55,0.3,0.35,0.6])

hold on
plot(traj(1,:),traj(2,:),'k*');
plot(traj2(1,1:152),traj2(2,1:152),'k*');

plot(g_grid.xs{1}(iv_RAS), g_grid.xs{2}(iv_RAS),'g.', 'MarkerSize',3)
visSetIm(g_p,goal_p,'b',0)
visSetIm(g_p,goal2_p,'k',0)
visSetIm(g_p,obs_p,'r',0)
visSetIm(g_p,clvf_p,'m',0)
visSetIm(g_p,ra(:,:,end),'k',0)


set(gca,'yTick',[-3:4:5]);
set(gca,'xTick',[-5:5:5]);
xlim([-5,5]);
ylim([-3,5]);
zx1 = get(gca,'ZTickLabel');
set(gca,'ZTickLabel',zx1,'fontsize',25);

lg1 = legend([OBS,GOAL,GOAL2,CLVF,TRAJ],...
    {'obstacle, ', 'target1, ', 'target2, ','$\mathcal I_m$, ','traj'}, ...
    'Interpreter', 'latex', 'FontSize', 18 , 'Orientation','horizontal' );