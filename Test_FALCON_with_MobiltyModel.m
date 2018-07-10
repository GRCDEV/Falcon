% SAMPLE SCRIPT FOR TESTING THE FALCON MODEL USING MOBILITY MODELS
%
% Developed by Enrique Hernandez-Orallo, Grupo de Redes de Computadores, Universitat Politecnica de Valencia, 2018.
%
% See paper: Enrique Hernandez-Orallo, Juan Carlos Cano, Carlos T. Calafate, Pietro Manzoni, 
% "FALCON: A New Approach for the Evaluation of Opportunistic Networks", Ad Hoc Networks 2018



% Type of mobility
RWM = 1;
RCM = 0;

PEDESTRIAN_MODEL = true;
VANET_MODEL = false;
SIMPLE_MODEL = false;


GRAPH_DYNAMIC = true;


%
% Simple model
%
if SIMPLE_MODEL
Sim_time = 4000;
dt = 10;  % 10s -> a time to move in the cell....
R = 10; C = 10;  % MxM cell.
CellSide = 10; % Cell side size (meters). In this case bluetooth
N = 3;  % Number of nodes
MOB_MODEL=RWM;  % Mobility model (RCM or RWM)
% Configuration for RCM
p_mov = 0.2; % Probability of movement in each step....
% Configuration for RWM
Speed_Interval = [0.2 2.2];  %(m/s)
Pause_Interval = [0 10];      %pause time (s)
Walk_Interval = [2.00 20.00]; %walk time (s)
Q = 30;
Cg =ones(R,C);
Cg = [ 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0;
       0 0 1 1 1 1 1 0 0 0;
       0 0 1 1 1 1 1 0 0 0;
       0 0 1 1 1 1 1 0 0 0;
       0 0 1 1 1 1 1 0 0 0;
       0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0;];
      
end


%
% Simple Bluetooth sample.
%
if PEDESTRIAN_MODEL
Sim_time = 8000;
dt = 10;  % 10s -> a time to move in the cell....
R = 15; C=20;  % RxC cell.
CellSide = 10; % Cell side size (meters). In this case bluetooth
N = 10;  % Number of nodes
MOB_MODEL=RWM;  % Mobility model (RCM or RWM)
% Configuration for RCM
p_mov = 0.2; % Probability of movement in each step....
% Configuration for RWM
Speed_Interval = [0.2 2.2];  %(m/s)
Pause_Interval = [0 10];      %pause time (s)
Walk_Interval = [2.00 20.00]; %walk time (s)
Q = R*C/4;
Cg =rand(R,C);
Cg =ones(R,C);
Cg = [ 0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0; 
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0;
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0; 
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0;
       0 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 0 0 1 0;
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ];

end



%
% VANET sample 10x10 km, 1000 vehicles, WiFi (100m), 12 hours.
%
if VANET_MODEL
    CellSide = 100; 
    R=100; C=100; N=1000;
    dt = 20;Sim_time = 3600;
    MOB_MODEL=RWM;
    Speed_Interval = [0.2 2.2];  %(m/s)
    Pause_Interval = [0 10];      %pause time (s)
    Walk_Interval = [20 200.00]; %walk time (s)
    Q=1000;
    Cg = ones(R,C);
end        

TotalArea = CellSide^2*R*C; % meters^2
item_x = ceil(Q/2)+1;   % Number of item to evaluate its diffusion


fprintf('Starting ... Cell Grid model\n');
fprintf('  Simulation: Max_T = %d, dt = %5.3f\n', Sim_time, dt); 
fprintf('  Scenario  : Rows= %d Cols = %d CellSide = %5.1f (meters), Area = %5.0f (m^2) NumItems = %6d\n', R, C , CellSide, TotalArea, Q); 
fprintf('  Nodes     : N= %d\n', N); 

t_start=tic;

% ONE: Generate items.
fprintf('Gen. Items... ');
Dp = FALCON_GenItemsRand(R,C,Q);
% Dp = [1, 1]; Q = 1;item_x = 1;

% SECOND: Generate nodes movements with respect the cell grid. 
fprintf('Gen. nodes movements... ');
if MOB_MODEL==RCM
    [Ng,Np]=FALCON_RCM_GenPoints(N,R,C,p_mov,Sim_time,dt); 
else
    RWM_input = struct('V_POSITION_X_INTERVAL',[0,C*CellSide],...%(m)
                 'V_POSITION_Y_INTERVAL',[0,R*CellSide],...%(m)
                 'V_SPEED_INTERVAL',Speed_Interval,...%(m/s)
                 'V_PAUSE_INTERVAL',Pause_Interval,...%pause time (s)
                 'V_WALK_INTERVAL',Walk_Interval,...%walk time (s)
                 'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
                 'SIMULATION_TIME',Sim_time,...%(s)
                 'NB_NODES',N);
    RWM_mobility = Generate_Mobility(RWM_input);    
    
    [Ng,Np]=FALCON_RWM_GenPoints(RWM_mobility,N,CellSide,Sim_time,dt);
end 
fprintf('OK!');

% THIRD: Simulate
fprintf(' Start sim... ');
% Initial state. No node has items.
X0 = zeros(N,Q); 
[Dn,Di,Ct,Ic,Im] = FALCON_Simulate(N,R,C,Np,Dp,Q,Cg,X0,Sim_time,dt,item_x, GRAPH_DYNAMIC);

fprintf('DONE %f (seconds)\n  ', toc(t_start));

TT = 0:dt:Sim_time;

RatioData = sum(Dn)/(N*Q);

figure;

subplot(2,2,1);
plot(TT,RatioData);
title('Ratio of data collected');
xlabel('time (s)');
ylabel('Ratio');

subplot(2,2,2);
plot(TT,Di);
title('Diffusion of center item');
xlabel('time (s)');
ylabel('Nodes');

subplot(2,2,4);
plot(TT(2:end),Ct(2:end));
title('Number of contacts');
xlabel('time (s)');
ylabel('Contacts');

subplot(2,2,3);
plot(TT,Ic,'r',TT,Im,'b');
title('Number of items');
legend('Items collected','Items interchanged');
xlabel('time (s)');
ylabel('Items');

% plot(TT,RatioData);
% xlabel('time (s)');
% ylabel('Ratio');
% set(gca,'fontsize',22);