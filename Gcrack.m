function crack_path_2 = Gcrack(Normal_stress, Number_of_Nodes, nodal_connectivity_values, sigma_t, Coordinates_Plate)

[rows,columns] = size(Normal_stress);

Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values);
% %Change to checking all nodes repittions to see if any one 
% %exceeds sigma_t
% 
% 
% 
% for i = 1:rows
%     for j = 1:columns
%         if Repitition_Remover(i,j) == 1
%             if abs(Normal_stress(i,j)) > sigma_t
%                 crack_path_1(i,j) = Normal_stress(i,j);
%             end
%         end
%     end
% end

crack_path_1 = Node_Highest_Stress_Identifier(nodal_connectivity_values,Normal_stress,sigma_t);

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


for k = 1:Number_of_Nodes %going over the nodes
    if crack_path_2(k,2) ~= 0 %chekcing if the point is overstressed
        Length = 0;
        for l = k-4:k+4 %Looping over nearby nodes
            if l > 1 && l < Number_of_Nodes
            if crack_path_2(l,2) ~= 0 %Checking if the nearby nodes are overstressed
                    % Finding Differences and adding them up
                    Length = Length + (0.5) * sqrt((Coordinates_Plate(k,1) - Coordinates_Plate(l,1))^2 +...
                            (Coordinates_Plate(k,2) - Coordinates_Plate(l,2))^2); 
            end
            end
        end
        crack_path_2(k,3) = Length; %%assisning additions to the matrix
    end
end



end