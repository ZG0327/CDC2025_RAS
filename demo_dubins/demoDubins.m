close all
clear 
clc


%% Grid
grid_min = [-5; -5; -pi]; % Lower corner of computation domain
grid_max = [5; 5; pi];    % Upper corner of computation domain
N = [71; 71; 51];         % Number of grid points per dimension
pdDims = 3;               % 3rd dimension is periodic
g = createGrid(grid_min, grid_max, N, pdDims);

%% time vector
t0 = 0;
tMax = 5;
dt = 0.1;
tau = t0:dt:tMax;

%% problem parameters
uMode = 'min';
dMode = 'max';



%% Dynamics
params.v = 1; % Velocity of the Dubins car
params.pxd = 0; % Desired postion
params.pyd = 0; % Desired velocity
params.u_max = pi; % maximum control input
params.u_min  = -pi; % minimum control input 

wRange = [ params.u_min , params.u_max ];
dRange = {[-0.2;-0.2;0] ; [0.2; 0.2; 0]};
speed = params.v;
dCar = DubinsCar([0, 0, 0], wRange, speed, dRange);

gamma = 0.2;

%% Pack problem parameteres
schemeData.grid = g;
schemeData.dynSys = dCar;
schemeData.accuracy = 'high'; %set accuracy
schemeData.uMode = uMode;
schemeData.dMode = dMode;
schemeData.clf.gamma = gamma;

%% Compute the CLVF using algorithm 1. 
% First find the smallest value for the value function with gamma = 0
% Test with different cost functions:
% 1: l(x) = ||x||_2, 
% 2: l(x) = ||x||_infty
% 3: l(x) = x'Qx
% 
% data01 = shapeRectangleByCorners(g, [-1;-1;-inf], [3;2;inf]);
% data02 = shapeCylinder(g, 3, [-2; -2; 0], 1); 
% data03 = shapeCylinder(g, 3, [3; -3; 0], 1); 
% 
% obs = min(data01,data02);
goal = shapeCylinder(g, 3, [3.5; 3.5; 0], 0)-0.8;

% data0 = importdata("fine_data\V0_2norm.mat");
% data0 = data0.value - 1.9;
% Q = [2,0,0;0,1,0;0,0,1];
% data0_1 = QuadCost(g,Q);

% goal = shapeCylinder(g, [], [0; 0], 0.5);
% data0 = max(goal,-obs);
data0 = goal;
%%
% data0: initial value function
% cost: cost function
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.valueFunction = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.deleteLastPlot = true; 
HJIextraArgs.convergeThreshold = 0.01;
HJIextraArgs.stopConverge = 1;
HJIextraArgs.stopDiverge = 1;
HJIextraArgs.divergeThreshold = 10;
HJIextraArgs.keepLast = 1;
HJIextraArgs.makeVideo = 0;
HJIextraArgs.targetFunction = goal;
% HJIextraArgs.obstacleFunction = obs;
HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; 
HJIextraArgs.visualize.plotData.projpt = {'min'}; 
HJIextraArgs.visualize.viewAngle = [70,25];
HJIextraArgs.visualize.xTitle = 'x';
HJIextraArgs.visualize.yTitle = 'y';
% HJIextraArgs.visualize.targetFunction = 1;
% HJIextraArgs.visualize.obstacleFunction = 1;
HJIextraArgs.visualize.plotColorTF = 'm';
% HJIextraArgs.visualize.plotColorOF = 'r';
% HJIextraArgs.visualize.plotAlphaOF = 0.5;
HJIextraArgs.visualize.plotAlphaTF = 0.5; 
HJIextraArgs.videoFilename = 'original_RA';

[V, tau1, ~] = ...
  HJIPDE_ZGsolve(data0, tau, schemeData, 'minCLF', HJIextraArgs);
% maxVwithL

