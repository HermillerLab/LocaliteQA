function new_rotmax = rem_NaN_rotmax(rotmax)
% This function removes NaN values from the last row of a rotation matrix
% and replaces them with 0

% Reshape the input rotation matrix
new_rotmax = reshape(rotmax, size(rotmax, 1), 16);

% Loop through each row of the reshaped matrix
for i = 1:size(new_rotmax, 1)
    % Check elements in columns 13 to 15 
    for j = 13:15
        % If an element is NaN, replace it with 0
        if isnan(new_rotmax(i, j))
            new_rotmax(i, j) = 0.0;
        end
    end
end

end