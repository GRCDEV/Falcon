function [Dp] = FALCON_GenItemsRand(R,C,num_items)
    % Generate the information grid with num_items randomly positioned.
    % Only one item by cell....
    % INPUT
    %    R : number of rows of the grid. 
    %    C : number of col of the grid.
    %    num_items : number of items to generate
    % OUTPUT
    %    Dp : Items position. A matrix of num_item x 2, with the position (y,x)
    %         of the items.
    Dg = zeros(R,C); 
    Dp = zeros(num_items,2);
    item = 1;
    if num_items == R*C
        for yy=1:R
            for xx=1:C
                Dp(item,:) = [yy,xx]; 
                item = item + 1;
            end
        end
    else
        % Generate the position randomly
        for ii=1:num_items
            x = randi(C); y = randi(R);
            while Dg(y,x) > 0
                x = randi(C); y = randi(R);              
            end
            Dg(y,x) = ii;
            Dp(ii,:) = [y,x];
        end
    end
end