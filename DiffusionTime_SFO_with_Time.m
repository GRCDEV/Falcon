% SAMPLE SCRIPT FOR TESTING THE FALCON MODEL USING THE SFO TRACE
% Developed by Enrique Hernandez-Orallo, Grupo de Redes de Computadores, Universitat Politecnica de Valencia, 2018.
%
% See paper: Enrique Hernandez-Orallo, Juan Carlos Cano, Carlos T. Calafate, Pietro Manzoni, 
% "FALCON: A New Approach for the Evaluation of Opportunistic Networks", Ad Hoc Networks 2018
clear;
load('SFO_trace_2008_5_18_area.mat');
GRAPH_DYNAMIC = false;

% Approximate conversion from degrees to meter for San Francisco Latitudes
LatDeg2m = 111319;
LonDeg2m = 88055;

Lon_min = -122.47;
Lon_max = -122.37;
Lat_min = 37.72;
Lat_max = 37.82;

XMax = (Lon_max - Lon_min)*LonDeg2m; 
YMax = (Lat_max - Lat_min)*LatDeg2m; 
CellSide = 140;
R = ceil(YMax/CellSide);
C = ceil(XMax/CellSide);

dt = 1;
c_radio = 100;

% From 7am to 11am 
start_time = datetime(2008,5,18,7,0,0);
posix_start_time = posixtime(start_time);
end_time = datetime(2008,5,18,12,0,0);
posix_end_time = posixtime(end_time);


trace = GPSTrace_extract_from_interval(SFO_trace_2008_5_18_area,posix_start_time,posix_end_time);

Sim_time = 3600*4;
Cg = ones(R,C);



TotalArea = CellSide^2*R*C; % meters^2
item_x = 1;   % Number of item to evaluate its diffusion


fprintf('Starting ... Cell Grid model\n');
fprintf('  Simulation: Max_T = %d, dt = %5.3f\n', Sim_time, dt); 
fprintf('  Scenario  : Rows= %d Cols = %d CellSide = %5.1f (meters), Area = %5.0f (m^2)\n', R, C , CellSide, TotalArea); 

t_start=tic;



% Generate items.
fprintf('Gen. Items... ');
% At fixed points in order to avoid the sea.
Dp = [  12 45; 27 45; 40 45; 53 45; 65 45;
       12 32; 27 32; 40 32; 53 32; 65 32;
       12 19; 27 19; 40 19; 53 19; 65 19; 
        12 7; 27  7; 40  7; 53  7; 65  8];
% El item 67 7 se comporta muy mal... lo muevo a las 65 8
Sp = size(Dp);
Q = Sp(1);
fprintf('Q = %d ', Q);

trace(:,2) = trace(:,2)-posix_start_time; 


% Convert to meters from Lat_min and Lon_min
trace(:,3) = (trace(:,3)-Lat_min)*LatDeg2m; 
trace(:,4) = (trace(:,4)-Lon_min)*LonDeg2m;


%%%%%%%%%%%%%%%%%%%%
%  FALCON Model
%%%%%%%%%%%%%%%%%%%%

% First: Generate nodes movements with respect the cell grid. 
fprintf('Gen. nodes movements... ');
   
[Ng,Np]=FALCON_GPS_GenPoints(trace,CellSide, Sim_time,dt);
N = length(Ng);
fprintf('OK!... Nodes %d\n', length(Ng));

TT = 0:dt:Sim_time;

T_setup = 10;


X0 = zeros(N,Q); 
Ratio = 0.99;

jj = 1;
for T_item = 0:1:2 % 15
    
    T_i(jj) = T_item;

    fprintf(' Start sim FALCON with time... T_item = %f ', T_item);
    % Initial state. No node has items.
    
    T_items = ones(1,Q)*T_item;
    
    T_setup = 0;

    [ DnC,DiC,Ct,Ic,Im] = FALCON_Simulate_with_Time(N,R,C,Np,Dp,Q,Cg,X0,T_setup,T_items, Sim_time,dt,item_x, GRAPH_DYNAMIC, false);
    RatioData = sum(DnC)/(N*Q);
    
    % Obtain diffusion time
    for i = 1: length(RatioData)
        if RatioData(i) >= Ratio
            break;
        end
    end    
    DT0(jj) = TT(i);
    
    T_setup = 10;

    [ DnC,DiC,Ct,Ic,Im] = FALCON_Simulate_with_Time(N,R,C,Np,Dp,Q,Cg,X0,T_setup,T_items, Sim_time,dt,item_x, GRAPH_DYNAMIC, false);
    RatioData = sum(DnC)/(N*Q);
    
    % Obtain diffusion time
    for i = 1: length(RatioData)
        if RatioData(i) >= Ratio
            break;
        end
    end    
    DT1(jj) = TT(i);
    
    T_setup = 20;

    [ DnC,DiC,Ct,Ic,Im] = FALCON_Simulate_with_Time(N,R,C,Np,Dp,Q,Cg,X0,T_setup,T_items, Sim_time,dt,item_x, GRAPH_DYNAMIC, false);
    RatioData = sum(DnC)/(N*Q);
    
    % Obtain diffusion time
    for i = 1: length(RatioData)
        if RatioData(i) >= Ratio
            break;
        end
    end    
    DT2(jj) = TT(i);
    
    T_setup = 30;

    [ DnC,DiC,Ct,Ic,Im] = FALCON_Simulate_with_Time(N,R,C,Np,Dp,Q,Cg,X0,T_setup,T_items, Sim_time,dt,item_x, GRAPH_DYNAMIC, false);
    RatioData = sum(DnC)/(N*Q);
    
    % Obtain diffusion time
    for i = 1: length(RatioData)
        if RatioData(i) >= Ratio
            break;
        end
    end    
    DT3(jj) = TT(i);



    jj = jj + 1;

    fprintf('DONE %f (seconds)\n', toc(t_start));

end



plot(T_i,DT0,'b-o');
hold on;
plot(T_i,DT1,'r--x');
hold on;
plot(T_i,DT2,'b-*');
hold on;
plot(T_i,DT3,'r--+');
hold on;

xlabel('Item transmission time (s)');
ylabel('Diffusion time (s)');
set(gca,'fontsize',22);
legend('\tau=0','\tau=10','\tau=20','\tau=30');
