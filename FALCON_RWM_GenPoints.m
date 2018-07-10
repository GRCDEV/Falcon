function [Ng,Np] = FALCON_RWM_GenPoints(RWM_mobility,N, CellSide, Sim_time, dt)
    % From a generated Randow Waypoint Model, generate the cell movements.
    % It also return the initial position of all nodes.
    % INPUT: (al values in m, s, and m/s.
    %    RWM_mobility: mobility structure generate from Generate_Mobility
    %    N : Number of nodes
    %    CellSide: side of cell in meters.
    %    Sim_time : simulation time (s)
    %    dt : time step used to generate the points in s.
    % RETURNS:
    %    Ng: Matrix with the initial position of the nodes (rows, cols)
    %    Np: an array of structures. Each element corresponds to the
    %    lists of positions for the nodes.
    %               Np(NODE).PosXY, where PosXY is a matrix with the
    %               position x and y for each time interval.    
    v_t = 0:dt:Sim_time;
    
    for nodeIndex = 1:N
        %Simple interpolation (linear) to get the position, anytime.
        %Remember that "interp1" is the matlab function to use in order to
        %get nodes' position at any continuous time.
        Np(nodeIndex).PosYX(1,:) = ceil(interp1(RWM_mobility.VS_NODE(nodeIndex).V_TIME,RWM_mobility.VS_NODE(nodeIndex).V_POSITION_Y,v_t)/CellSide);
        Np(nodeIndex).PosYX(2,:) = ceil(interp1(RWM_mobility.VS_NODE(nodeIndex).V_TIME,RWM_mobility.VS_NODE(nodeIndex).V_POSITION_X,v_t)/CellSide);

    end
    Ng = zeros(N,2);
    for ii = 1:size(Ng)
        y = Np(ii).PosYX(1,1);
        x = Np(ii).PosYX(2,1);  
        Ng(ii,1) = y; Ng(ii,2) = x; 
    end     

end