function nodes_configuration = Node_Configurator(Coordinates)

[rows,columns] = size(Coordinates);
nodes_configuration = 0;

Row_Check = 0;
current_Row = 1;
row_index = 2;
    for i = 2:rows
       if Coordinates(i,1) ~= Coordinates(i-1,1) && Row_Check == 0
           nodes_configuration = current_Row;
           current_Row = 0;
           Row_Check = 1;
           row_index = 1;
       else if Coordinates(i,1) ~= Coordinates(i-1,1)
           nodes_configuration = [nodes_configuration;current_Row];
           row_index = 1;
       end
       end
       current_Row(row_index) = i;
       row_index = row_index + 1;
    end

nodes_configuration = [nodes_configuration;current_Row];

end