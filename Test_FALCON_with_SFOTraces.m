% SAMPLE SCRIPT FOR TESTING THE FALCON MODEL USING THE SFO TRACE
%
% Developed by Enrique Hernandez-Orallo, Grupo de Redes de Computadores, Universitat Politecnica de Valencia, 2018.
%
% See paper: Enrique Hernandez-Orallo, Juan Carlos Cano, Carlos T. Calafate, Pietro Manzoni, 
% "FALCON: A New Approach for the Evaluation of Opportunistic Networks", Ad Hoc Networks 2018

load('SFO_trace_2008_5_18_area.mat');
GRAPH_DYNAMIC = true;

% Approximate conversion from degrees to meter for San Francisco Latitudes
LatDeg2m = 111319;
LonDeg2m = 88055;


CellSide = 130;
R=86; C=68; 
dt = 1;

% From 7am to 11am 
start_time = datetime(2008,5,18,7,0,0);
posix_start_time = posixtime(start_time);
end_time = datetime(2008,5,18,12,0,0);
posix_end_time = posixtime(end_time);


trace = GPSTrace_extract_from_interval(SFO_trace_2008_5_18_area,posix_start_time,posix_end_time);

Sim_time = 3600*3;
Cg = ones(R,C);

TotalArea = CellSide^2*R*C; % meters^2
item_x = 1;   % Number of item to evaluate its diffusion


fprintf('Starting ... Cell Grid model\n');
fprintf('  Simulation: Max_T = %d, dt = %5.3f\n', Sim_time, dt); 
fprintf('  Scenario  : Rows= %d Cols = %d CellSide = %5.1f (meters), Area = %5.0f (m^2)\n', R, C , CellSide, TotalArea); 

t_start=tic;

% ONE: Generate items.
fprintf('Gen. Items... ');
% Dp = CGM_GenItemsRand(R,C,Q);
% At fixed points in order to avoid the sea.
Dp = [ 30 30; 40 30; 50 30; 60 30; 70 30;
       30 40; 40 40; 50 40; 60 40; 70 40;
       30 20; 40 20; 50 20; 60 20; 70 20
       10 2; 7 2; 30 2; 40 2; 50 2];
Q = length(Dp);
fprintf('Q = %d ', Q);

% SECOND: Generate nodes movements with respect the cell grid. 
fprintf('Gen. nodes movements... ');

trace(:,2) = trace(:,2)-posix_start_time; 

Lon_min = -122.47;
Lon_max = -122.37;
Lat_min = 37.72;
Lat_max = 37.82;

% Convert to meters from Lat_min and Lon_min
trace(:,3) = (trace(:,3)-Lat_min)*LatDeg2m; 
trace(:,4) = (trace(:,4)-Lon_min)*LonDeg2m; 


[Ng,Np]=FALCON_GPS_GenPoints(trace,CellSide, Sim_time,dt);
N = length(Ng);
fprintf('OK!... Nodes %d\n', length(Ng));


% THIRD: Simulate
fprintf(' Start sim... ');
% Initial state. No node has items.
X0 = zeros(N,Q); 
[Dn,Di,Ct,Ic,Im] = FALCON_Simulate(N,R,C,Np,Dp,Q,Cg,X0,Sim_time,dt,item_x, GRAPH_DYNAMIC, false);

fprintf('DONE %f (seconds)\n  ', toc(t_start));

TT = 0:dt:Sim_time;

RatioData = sum(Dn)/(N*Q);

subplot(2,2,1);
plot(TT,RatioData);
title('Ratio of data collected');
xlabel('time (s)');
ylabel('Ratio');

subplot(2,2,2);
plot(TT,Di);
title('Diffusion of center item');
xlabel('time (s)');

subplot(2,2,4);
plot(TT(2:end),Ct(2:end));
title('Number of contacts');
xlabel('time (s)');

subplot(2,2,3);
plot(TT,Ic,'r',TT,Im,'b');
title('Number of items');
legend('Items collected','Items interchanged');
xlabel('time (s)');

% plot(TT,RatioData);
% xlabel('time (s)');
% ylabel('Ratio');
% set(gca,'fontsize',22);