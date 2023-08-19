function nodes_configuration = Node_Configurator(Coordinates)
%% Going over the coordinates file and creating a matrix with the node network so it can be used later to find the crack path


[rows,columns] = size(Coordinates);
nodes_configuration = 0;

%% Checking to see which column is changing
Column_Changing = 1;

if Coordinates(1,1) ~= Coordinates(2,1) %if the first column is being varied, 
    Column_Changing = 2; %It means the second column is the pivot
else
    Column_Changing = 1; %Or otherwise
end


Row_Check = 0;
current_Row = 1;
row_index = 2;
    for i = 2:rows
       if Coordinates(i,Column_Changing) ~= Coordinates(i-1,Column_Changing) && Row_Check == 0
           nodes_configuration = current_Row;
           current_Row = 0;
           Row_Check = 1;
           row_index = 1;
       else if Coordinates(i,Column_Changing) ~= Coordinates(i-1,Column_Changing)
           nodes_configuration = [nodes_configuration;current_Row];
           row_index = 1;
       end
       end
       current_Row(row_index) = i;
       row_index = row_index + 1;
    end

nodes_configuration = [nodes_configuration;current_Row];

nodes_configuration = flip(nodes_configuration);


end