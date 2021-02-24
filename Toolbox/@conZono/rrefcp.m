function [out] = rrefcp(obj)

A_rref = [obj.A obj.b]; %concatenates A and b matrices to compute rref

[m,n] = size(A_rref); % computes the dimension of A_rref.


[num, den] = rat(A_rref); % Checks if elements of A_rref are ratios of small integers.
rats = isequal(A_rref, num./den);

ind = [1:n]; % Stores the column index

for j = 1:m % Iterates through each of the rows.
    
    temp_max = abs(A_rref(j,j)); % Sets A_rref(j,j) as the element with maximum absolute value.
    imax = j;
    jmax = j;
    
    for k = j:m
        for l = j:n-1
            if (temp_max < abs(A_rref(k,l))) % Determines the element with maximum absolute value in the matrix
                % between rows (j - m) and columns (j - (n-1)) matrix
                % (Pivot element) (Leaves out b vector which does not need reduction)
                temp_max = abs(A_rref(k,l));
                imax = k; % Stores the row index of the element with maximum absolute value
                jmax = l; % Stores the corresponding column index
            end
        end
    end
    
    temp_jmax = A_rref(:,jmax); % Swaps the jmax^{th} column and j^th column
    A_rref(:,jmax) = A_rref(:,j);
    A_rref(:,j) = temp_jmax;
    
    temp_col_index = ind(jmax); % Stores column indices after column swap
    ind(jmax) = ind(j);
    ind(j) = temp_col_index;
    
    temp_imax = A_rref(imax,:); % Swaps the imax^{th} row and j^th row
    A_rref(imax,:) = A_rref(j,:);
    A_rref(j,:) = temp_imax;
    
    if(j ~= m)    % Checks if the last row is reached
        for i = j+1:m % Traverses through each of the rows from (j+1) to (m)
            A_rref(i,:) = A_rref(i,:) - (A_rref(i,j)/A_rref(j,j)).*A_rref(j,:); % Performs row reduction
        end
    end
end

for i =1:m % Traverses through each of the rows from 1 to m
    if (A_rref(i,i) ~= 0) % Checks for non-zero.
        A_rref(i,:) = A_rref(i,:)/A_rref(i,i); % Normalizes the matrix to get 1 in A_rref(i,i) element
    end
end

if rats
    [num, den] = rat(A_rref); % Returns rational numbers if true.
    A_rref = num./den;
end

obj0 = conZono;
obj0.A = A_rref(:,1:size(obj.A,2));
obj0.b = A_rref(:,size(obj.A,2)+1);
obj0.G = obj.G(:,ind(1:end-1));
obj0.c = obj.c;

out = obj0;
end