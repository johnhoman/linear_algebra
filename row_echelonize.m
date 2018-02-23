function S = row_echelonize(A)
    tic;
    if any(size(A) <= 1)
        S.time = toc;
        S.A = A;
        return;
    end

    [rows, cols] = size(A);
    
    row = 1;
    col = 1;
    while row < rows && col < cols
        pivot_row = find(A(:, col)' ~= 0 & (1:rows) >= row);
        if isempty(pivot_row)
            col = col + 1;
            continue;
        end
        pivot_row = pivot_row(1);
        
        if pivot_row > row
            A([row pivot_row], :) = A([pivot_row row], :);
        end
        
        pivot = A(row, col);
        for k = row + 1:rows
            front = A(k, col);
            A(k, :) = A(k, :) - front/pivot*A(row, :);
        end
        
        row = row + 1;
        col = col + 1;
    end
    S.time = toc;
    S.A = A;
end