function [Crakcs_All,crack_path_2] = GcrackPath(Normal_stress, Number_of_Nodes, nodal_connectivity_values, sigma_t, Coordinates_Plate,Purpose)

%% Function used to calculate the nodal Stresses from a matrix with element stresses
[rows,columns] = size(Normal_stress);

Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values);

crack_path_1 = Node_Highest_Stress_Identifier(nodal_connectivity_values,Normal_stress,sigma_t,Purpose);

crack_path_2 = zeros(21,3);

%% Getting Element Numbers too

for k = 1:Number_of_Nodes
    for i = 1:rows
        for j = 1:columns
            if nodal_connectivity_values(i,j) == k
                if Repitition_Remover(i,j) == 1
                    crack_path_2(k,1) = nodal_connectivity_values(i,j);
                    crack_path_2(k,2) = crack_path_1(i,j);
                end
            end
        end
    end
end

Cracks_All = {};
node_configuration = Node_Configurator(Coordinates_Plate);
node_configuration = node_configuration';
node_configuration = flipud(node_configuration);
node_configuration = fliplr(node_configuration);

for k = 1:Number_of_Nodes %going over the nodes
    if crack_path_2(k,2) ~= 0 %chekcing if the point is overstressed
        Crack_Path_Draw = [];
        Length = 0;
        [row_node_configuration,column_node_configuration] = size(node_configuration);
        [row,column] = find(node_configuration == k);
        %Left Right Initiation
        left = column - 1;
        right = column + 1;
        no_more_right = false;
        no_more_left = false;

        %Up Down Initiation
        above = row - 1;
        below = row + 1;
        no_more_above = false;
        no_more_below = false;


        Directions = {};

        Directions(1) = {[below,left]};
        Directions(2) = {[row,left]};
        Directions(3) = {[above,left]};

        Directions(4) = {[above,column]};
        Directions(5) = {[below,column]};

        Directions(6) = {[above,right]};
        Directions(7) = {[row,right]};
        Directions(8) = {[below,right]};

        for l = 1:8
            if Directions{l}(1)> 0 && Directions{l}(1) < row_node_configuration + 1 &&  Directions{l}(2)> 0 && Directions{l}(2) < column_node_configuration + 1

                Coordinates_Original = [Coordinates_Plate(node_configuration(row,column),1),Coordinates_Plate(node_configuration(row,column),2)];

                Node_at_check = node_configuration(Directions{l}(1),Directions{l}(2));
                Coordinates_Next = [Coordinates_Plate(Node_at_check,1),Coordinates_Plate(Node_at_check,2)];

                Length = Length + (0.5) * sqrt(((Coordinates_Original(k,1) - Coordinates_Next(1))^2 + ...
                    Coordinates_Original(k,2) - Coordinates_Next(2))^2);

                Crack_Path_Temp = [Node_at_check, Coordinates_Next(1), Coordinates_Next(2)];
                Crack_Path_Draw = [Crack_Path_Draw; Crack_Path_Temp];

            end
        end
        crack_path_2(k,3) = Length; %%assisning additions to the matrix
        [row_Cracks,column_Cracks] = size(Cracks_All);
        Cracks_All(row_Cracks + 1,1) = {Crack_Path_Draw};
    end
end



end