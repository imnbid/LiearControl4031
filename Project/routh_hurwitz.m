function routh_table = routh_hurwitz(coeffs)
    n = length(coeffs); % Order of polynomial
    routh_table = zeros(n, ceil(n/2)); % Initialize Routh array

    % Fill the first two rows
    routh_table(1, :) = coeffs(1:2:end); % Even index coefficients
    routh_table(2, 1:length(coeffs(2:2:end))) = coeffs(2:2:end); % Odd index

    % Compute remaining rows
    for i = 3:n
        for j = 1:size(routh_table, 2)-1
            num = routh_table(i-1, 1) * routh_table(i-2, j+1) - ...
                  routh_table(i-2, 1) * routh_table(i-1, j+1);
            den = routh_table(i-1, 1);
            if den == 0
                den = 1e-6; % Avoid division by zero
            end
            routh_table(i, j) = num / den;
        end

        % If the entire row is zero, replace with auxiliary equation
        if all(routh_table(i, :) == 0)
            order = n - i;
            aux_poly = [order:-2:0] .* routh_table(i-1, 1:order/2+1);
            routh_table(i, :) = [aux_poly, zeros(1, size(routh_table, 2) - length(aux_poly))];
        end
    end
end

