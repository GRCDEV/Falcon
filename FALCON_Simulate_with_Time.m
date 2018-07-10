function [Dn,Di,Ct,Ic,Im,DMI]= FALCON_Simulate_with_Time(N,R,C,Np,Dp,Q,Cg,X0,T_setup,T_items,Sim_time,dt,item_x,GRAPH_DYNAMIC,GRAPH_REVERSE)
% Evaluates the FALCON model with communications set-up time (this version
% is an extension of the FALCON_Simulate function considering this time.)
% Developed by Enrique Hernandez-Orallo, Grupo de Redes de Computadores, Universitat Politecnica de Valencia, 2018.
%
% See paper: Enrique Hernandez-Orallo, Juan Carlos Cano, Carlos T. Calafate, Pietro Manzoni, 
% "FALCON: A New Approach for the Evaluation of Opportunistic Networks", Ad Hoc Networks 2018
% 
% INPUT: 
%    N : number of nodes
%    R : number of rows of the grid.
%    C : number of columns of the grid.
%    Np: an array of structures. Each element corresponds to the
%    lists of positions for the nodes.
%    Dp : Items position. A matrix of num_item x 2, with the position (y,x)
%         of the items. Affects the interchanging of messages.
%    Q : Number of items
%    Cg : Communication grid. The probability of communication in a cell. 
%            From 0 (no communication) to 1 (full communication)
%    X0 : Initial state. Each i row corresponds to the items that have node i.
%    T_setup: connection setup time (the time required to start a
%    transmission). If T_setup < 0, then no time is considered in the interchange.
%    T_items: a vector of size N with the time required to transmit each item
%    Sim_time : simulation time (s)
%    dt : time step used to generate the points in the results.
%    item_x : numnber of item for evaluating its diffusion.
%    GRAPH_DYNAMIC: If true, displays a figure with the real-time
%    status of the diffusion,
%    REVERSE: If true, displays the grid without reverse (by default is
%    true)
% RETURNS:
%    Dn: Number of items in each node depending on time.
%         Is an array of Nx(Sim_time/dt) 
%    Di : Diffusion of item_x, that is the number of nodes that have it
%         depending on time
%    Ct : Number of contacts depending on time
%    Ic : Items collected.
%    Im : Items interchanged.
%    DMI: Density of messages interchanged by cell

% THIS FUNCION IS THE SAME THAT THE ONE IN CONTACT SIMULATION 
function [t_Trans, nm] = InterchangeItems(node_i,node_j,t_Contact,t)    
% This local function interchange items between node i and j with a
% duration of t_contact at time t.
% Interchange messages depending on time
% Returns 
%    t_Trans = time of transmission (or -1 if no transmission)
%    nm = number of messages interchanged.
% 
%  ***** Global variables used ***********************************
%
%      X: a matrix of NxQ size indicating if a node have an item.
%      node_trans: a vector of size N indicating if a node is transmitting
%      N, Q: number of nodes and items.
%      T_setup: connection setup time
%      T_items: a vector of size N with the time required to transmit each item
%
%    ---> Perfomance issues: to avoid passing it each time the function is called 
%  ***************************************************************

    nm = 0;
    t_Left = t_Contact; % Remaining time of the contact
    if T_setup < t_Left 
        t_Left = t_Left - T_setup;
        for q = 1:Q
            % If one of the nodes has the item.
            if (X(node_i,q) + X(node_j,q)) == 1
                if T_items(q) <= t_Left
                    X(node_i,q) = 1;
                    X(node_j,q) = 1;
                    nm = nm + 1;
                end
                t_Left = t_Left - T_items(q);
                if t_Left <= 0
                    break;  % No more items can be interchanged!!!
                end
            end
        end     
    else
        t_Left = 0; % Not enough time for the setup, but the time if consumed...
    end
    t_Trans = t_Contact - t_Left; % Final transmission time 
    if t_Trans > 0
        node_trans(node_i) = t_Trans+t;
        node_trans(node_j) = t_Trans+t;
    end
end

function P = GetPositionGrid(R,C)
    % Obtains a matrix with the number of nodes in each position of the grid.
    P = zeros(R,C);
    for i = 1:size(Ng)
        y = Ng(i,1); x = Ng(i,2);
        P(y,x) = P(y,x) + 1;
    end 
end

function [Dg] = GenItemsMatrix(R,C,Dp)
    % Generate the items grid. Items grid. A matrix of size RxC with the location of the items.
    %           y corresponds to the row, and x to the column.
    %         If y,x is > 0, then there is data to collect. 
    % Only used for optimizing the simulation WHEN max 1 item per cell.
    Dg = zeros(R,C);
    num_items = size(Dp);
    for ii=1:num_items
        Dg(Dp(ii,1),Dp(ii,2)) = ii;
    end
end


function lNg = MoveNodes(tt)
    % Move nodes to the next position using the points.
    % ****UPDATES Ng and node_exit ******
    % Return a list of all the nodes that change position from the last iteration
    % It also computes the time when it leaves the cell
    lNg = [];
    for ii = 1:N
        y = Np(ii).PosYX(1,tt);
        x = Np(ii).PosYX(2,tt);  
        
        if y ~= Ng(ii,1) || x ~= Ng(ii,2)
            Ng(ii,1) = y; Ng(ii,2) = x; 
            % Computes the exit time from the cell
            for tt2 = tt:TamTT
                if y ~= Np(ii).PosYX(1,tt2) || x ~= Np(ii).PosYX(2,tt2)
                    break;
                end
            end
            lNg = [lNg; [ii,y,x]]; 
            node_exit(ii) = TT(tt2);
        end
    end 
end

function lNg = EndCommNodes(t)
    % Check the nodes that have ended their communication up to time t
    % ****UPDATES node_trans ******
    % Return a list of all the nodes that ended their communication
    lNg = [];
    for ii = 1:N
        if node_trans(ii) > 0 && node_trans(ii) <= t
            node_trans(ii) = 0;
            y = Ng(ii,1);
            x = Ng(ii,2);  
            lNg = [lNg; [ii,y,x]];
        end
    end 
end


function UpdateData_fromLocation(lNg)
    % Update the data from the cell.
    % lNg = list of nodes (y,x) position, that moved since the last update
    % ***UPDATE*** the matrix D   
    for ii = 1:size(lNg)
        n = lNg(ii,1); y = lNg(ii,2); x = lNg(ii,3);
        if Dg(y,x) > 0
            X(n,Dg(y,x)) = 1;
        end
    end
end

function [NumContacts] = UpdateData_fromContacts(lNg)
    % Update the data from the contact of the cell.
    % lNg = list of nodes (y,x) position, that moved since the last update
    % ***UPDATE*** the matrix D  
    % Return the number of contacts.
    NumContacts = 0;
    for ii = 1:size(lNg)
        nn = lNg(ii,1); y = lNg(ii,2); x = lNg(ii,3);
        p_c = Cg(y,x);
        for jj = 1:N
            % Find  nodes in the same cell 
            if y == Ng(jj,1) && x == Ng(jj,2)
                if rand <= p_c
                    if node_trans(jj) == 0 && nn ~= jj
                        % Obtain the contact time and then performs the contact
                        t_Contact = min(node_exit(ii), node_exit(jj))-TT(tt); 
                        NumContacts = NumContacts + 1;
                        [ T_trans, nm] = InterchangeItems(nn,jj,t_Contact,TT(tt));
                        DMI(y,x) = DMI(y,x) + nm;  % For statistics
                        if  T_trans > 0
                            break;  %% If 0 it is assumed 0 transmition time so it can try with further contacts.
                        end
                    end
                end
            end
        end
    end
end

function Ds = Get_Ds(R,C,N,Ds)
    % Get data diffusion grid
    Ds = zeros(R,C);
    NumItems = size(Dp);
    for nn = 1:N   
        for ii = 1:NumItems
            Ds(Dp(ii,1),Dp(ii,2)) = Ds(Dp(ii,1),Dp(ii,2)) +X(nn,ii);
        end
    end 
end

function In = Get_ItemsNodes(N)
    % Obtain the number of items that have each node
    In = zeros(1,N);
    for nn = 1:N
        In(nn) = sum(X(nn,:));
    end
end



function GraphIt(t,Pg,Ds,In,NumItems) 

    persistent Xrange;
    persistent Yrange;
    persistent map_nodes; 
    persistent map_items;
    persistent map_data;
    persistent NN;
    persistent Aux_g;
    FONTSIZE = 18;
    N_max_per_pos = 5;
    if t==0
        Xrange = 1:(C+1);
        Yrange = 1:(R+1);
        NN = 1:(length(In)+1);
        map_nodes = flipud(bone(N_max_per_pos));
        map_items = flipud(pink(N));
        map_data = flipud(gray(NumItems+1));
        Aux_g = zeros(R+1,C+1); % We need to add an extra row and col to display the MxM grid....
    end 
    
    ax = subplot(4,4,[1,2,5,6,9,10]);
    sTitle = sprintf('Nodes movement t=%5.2f s', t);
    Aux_g(1:R,1:C) = Pg;
    pcolor(Xrange,Yrange,Aux_g); 
    colormap(ax, map_nodes);
    caxis([0 N_max_per_pos]);
    title(sTitle);   
    set(gca,'XAxisLocation','top');
    if GRAPH_REVERSE
        set(gca,'Ydir','reverse');
    end
    set(gca,'fontsize',FONTSIZE);

    ax = subplot(4,4,[3,4,7,8,11,12]);
    
    sTitle = sprintf('Items received t=%5.2f s', t);
    Aux_g(1:R,1:C) = Ds;
    pcolor(Xrange,Yrange,Aux_g); 
    colormap(ax, map_items);
    colorbar;
    caxis([0 N]);
    title(sTitle);  
    set(gca,'XAxisLocation','top');
    if GRAPH_REVERSE
        set(gca,'Ydir','reverse');
    end
    set(gca,'fontsize',FONTSIZE);

    
    ax = subplot(4,4,[13,14]);

    sVector = sprintf('%3d ', In);
    if length(sVector) >= 20
        sTitle = sprintf('Items per node');
    else
        sTitle = sprintf('Items per node [ %s ]', sVector);
    end
    pcolor(NN,1:2,[In,0;In,0]); 
    colormap(ax,map_data);
    caxis([0 NumItems]);
    title(sTitle);   
    set(gca,'YTick',[])
    set(gca,'YColor','w');
    set(gca,'fontsize',FONTSIZE-2);
    
    ax = subplot(4,4,15);
    
    sTitle = sprintf('Connection grid');
    Aux_g(1:R,1:C) = Cg;
    pcolor(Xrange,Yrange,Aux_g); 
    colormap(ax,map_data);
    caxis([0 1]);
    colorbar;
    title(sTitle);   
    set(gca,'YTick',[])
    set(gca,'YColor','w');
    set(gca,'XAxisLocation','top');
    if GRAPH_REVERSE
        set(gca,'Ydir','reverse');
    end
    set(gca,'fontsize',FONTSIZE-2);
    
    ax = subplot(4,4,16);
    
    sTitle = sprintf('Items Location #Items=%d', NumItems);
    Aux_g(1:R,1:C) = Dg;
    pcolor(Xrange,Yrange,sign(Aux_g)); 
    colormap(ax,map_data);
    caxis([0 1]);
    title(sTitle);   
    set(gca,'YTick',[])
    set(gca,'YColor','w');
    set(gca,'XAxisLocation','top');
    if GRAPH_REVERSE
        set(gca,'Ydir','reverse');
    end
    set(gca,'fontsize',FONTSIZE-2);

end


% Status at each period. Each i row corresponds to the items that have node i.
X = X0; 
% Ng: Matrix with the current position of the nodes (rows, cols)
%   Initially is set to zero. Is updated in MoveNodes
Ng = zeros(N,2);

% End time of current transmission (0 if no transmission)...
node_trans = zeros(1,N);
% End time of current cell for all nodes  
node_exit = zeros(1,N);

% Density of interchanged messages
DMI = zeros(R,C);

Dg = GenItemsMatrix(R,C,Dp);  % Only in the case of max 1 items by cell


TT = 0:dt:Sim_time;

Dn = zeros(N,length(TT));
Di = zeros(1,length(TT));
Ct = zeros(1,length(TT));
Ic = zeros(1,length(TT));
Im = zeros(1,length(TT));

if nargin == 14 
    GRAPH_REVERSE = true;
end

Dtotal_ant = 0;
fprintf('       ');
TamTT = length(TT);
for tt = 1:TamTT
    fprintf('\b\b\b\b\b\b\b%6.2f%%',100*tt/TamTT);
    
    % First: obtain the nodes that change of cell
    lNg = MoveNodes(tt);
    
    % Second: check if these nodes can pick a new item....
    UpdateData_fromLocation(lNg);
    Ic(tt) = sum(sum(X))-Dtotal_ant; % New items collected
    
    % Third : obtain the nodes that ended their transmission
    lNg2 = EndCommNodes(TT(tt));
    
    % Fourth: All these node are the ones that must be checked for Contacts
    lNg = [ lNg; lNg2]; 
    Ct(tt) = UpdateData_fromContacts(lNg);
    
    % Finally obtain statistics
    Dn(:,tt) = Get_ItemsNodes(N);
    Dtotal = sum(Dn(:,tt)); 
    Im(tt) = Dtotal - Ic(tt) - Dtotal_ant; 
    Di(tt) = sum(X(:,item_x));
    Dtotal_ant = Dtotal;
    
    % Plot is optional
    if GRAPH_DYNAMIC
        Pg = GetPositionGrid(R,C);
        Ds = Get_Ds(R,C,N,Dp);
        In = Get_ItemsNodes(N);
        GraphIt(TT(tt),Pg,Ds,In,Q);
        pause(0.1);
    end 
    
end
fprintf('\b\b\b\b\b\b\b');





end
