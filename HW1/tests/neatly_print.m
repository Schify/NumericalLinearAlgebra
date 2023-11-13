function neatly_print(xs)
    [num_rows, num_cols] = size(xs);

    for i = 1:num_cols
        x = xs(:, i);
        for j = 1:num_rows
            fprintf("%.7s ", x(j))
            if j < num_rows
                fprintf("& ")
            end
        end
        fprintf("\n");
    end
end
