function [Ng,Np] = FALCON_GPS_GenPoints(GPS_trace, CellSide, Sim_time, dt)
    % From a GPS trace, generate the cell movements.
    % It also return the initial position of all nodes.
    % INPUT: (al values in m, s, and m/s.)
    %    GPS_trace: Trace is in the following format :
    %        Taxi ID, Time in UNIX epoch format, Latitude (m), Longitude (m), Occupied
    %    CellSide: side of cell in meters.
    %    Lat_min:  Latitude base
    %    Lon_min:  Longitude base
    %    Sim_time : simulation time (s)
    %    dt : time step used to generate the points in s.
    % RETURNS:
    %    Ng: Matrix with the initial position of the nodes (row, col)
    %    Np: an array of structures. Each element corresponds to the
    %    lists of positions for the nodes.
    %               Np(NODE).PosXY, where PosXY is a matrix with the
    %               position x and y for each time interval.  
    v_t = 0:dt:Sim_time;
    
    GPS_trace = sortrows(GPS_trace,1);   % Sort by node number

    nodeIndex = 0;
    i_start = 1;
    taxiID = GPS_trace(1,1);
    num_points = length(GPS_trace); 
    for i = 1:num_points
        if taxiID ~= GPS_trace(i,1) || i == num_points
            nodeIndex = nodeIndex + 1;
            V_TIME = GPS_trace(i_start:i-1,2);
            V_POSITION_Y = GPS_trace(i_start:i-1,3);
            V_POSITION_X = GPS_trace(i_start:i-1,4);
            % Assure that first and last time has position
            if V_TIME(1) > 0
                V_TIME = [0; V_TIME; ]; 
                V_POSITION_Y = [V_POSITION_Y(1); V_POSITION_Y];
                V_POSITION_X = [V_POSITION_X(1); V_POSITION_X];
            end
            
            if V_TIME(end) < Sim_time
                V_TIME = [V_TIME; Sim_time]; 
                V_POSITION_Y = [V_POSITION_Y; V_POSITION_Y(end)];
                V_POSITION_X = [V_POSITION_X; V_POSITION_X(end)];
            end
           
            
            %Simple interpolation (linear) to get the position, anytime.
            %Remember that "interp1" is the matlab function to use in order to
            %get nodes' position at any continuous time.
            Np(nodeIndex).PosYX(1,:) = ceil(interp1(V_TIME,V_POSITION_Y,v_t)/CellSide+0.001);
            Np(nodeIndex).PosYX(2,:) = ceil(interp1(V_TIME,V_POSITION_X,v_t)/CellSide+0.001);
          
            taxiID = GPS_trace(i,1);
            i_start = i;
            
        end

    end
    Ng = zeros(nodeIndex,2);
    for ii = 1:nodeIndex
        y = Np(ii).PosYX(1,1);
        x = Np(ii).PosYX(2,1);  
        Ng(ii,1) = y; Ng(ii,2) = x; 
    end     

end