function [out, out_ind] = rrefcp(obj)
obj0 = copy(obj);
obj0.getDimensions;

[num, den] = rat(obj0.A); % Checks if elements of A are ratios of small integers.
rats = isequal(obj0.A, num./den);

out_ind = [1:obj0.nG]; % Stores the column index

    for j = 1:obj0.nC % Iterates through each of the rows.
    
        temp_max = abs(obj0.A(j,j)); % Sets A(j,j) as the element with maximum absolute value.
        imax = j;
        jmax = j;
        
         for k = j:obj0.nC
            for l = j:obj0.nG-1
                if (temp_max < abs(obj0.A(k,l))) % Determines the element with maximum absolute value in the matrix 
                % between rows (j - m) and columns (j - (n-1)) matrix
                % (Pivot element) (Leaves out b vector which does not need reduction)
                temp_max = abs(obj0.A(k,l));
                imax = k; % Stores the row index of the element with maximum absolute value
                jmax = l; % Stores the corresponding column index
                end
            end
         end

    temp_jmax = obj0.A(:,jmax); % Swaps the jmax^{th} column and j^th column 
    obj0.A(:,jmax) = obj0.A(:,j);
    obj0.A(:,j) = temp_jmax;    
    
    temp_col_index = out_ind(jmax); % Stores column indices after column swap
    out_ind(jmax) = out_ind(j);
    out_ind(j) = temp_col_index;
    
    temp_imax = obj0.A(imax,:); % Swaps the imax^{th} row and j^th row
    obj0.A(imax,:) = obj0.A(j,:);
    obj0.A(j,:) = temp_imax;
        
        if(j ~= obj0.nC)    % Checks if the last row is reached
            for i = j+1:obj0.nC % Traverses through each of the rows from (j+1) to (m) 
                obj0.A(i,:) = obj0.A(i,:) - (obj0.A(i,j)/obj0.A(j,j)).*obj0.A(j,:); % Performs row reduction
            end 
        end
    end

    for i =1:obj0.nC % Traverses through each of the rows from 1 to m
        if (obj0.A(i,i) ~= 0) % Checks for non-zero.
        obj0.A(i,:) = obj0.A(i,:)/obj0.A(i,i); % Normalizes the matrix to get 1 in A(i,i) element
        end
    end

    if rats
        [num, den] = rat(obj0.A); % Returns rational numbers if true.
        obj0.A = num./den;
    end
	out = obj0;
end