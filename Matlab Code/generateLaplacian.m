function [ T ] = generateLaplacian( I, T_est)
%GENERATELAPLACIAN2 Summary of this function goes here
%   Detailed explanation goes here

    dimI = size(I);

    % Taking a box around the pixel
    win_size = 1;
    win_pixels = 9;
    
    % As per equation in paper when computing the laplacian
    % U is to be added to the window covariance
    U = .000001 ./win_pixels.*eye(3);
    
    windowIndicies = 1:dimI(1)*dimI(2);
    windowIndicies = reshape(windowIndicies,dimI(1),dimI(2));
    
    totalElements = win_pixels^2 * ( dimI(1) - 2 ) * ( dimI(2) - 2 );
    
    indicies_x = ones(1,totalElements);
    indicies_y = ones(1,totalElements);
    elements = zeros(1,totalElements);
    
    count = 0;
    for i = (1+win_size):(dimI(2)-win_size)
        for j = (1+win_size):(dimI(1)-win_size)
            
            % Get the window around i and j
            rangeI = i-win_size:i+win_size;
            rangeJ = j-win_size:j+win_size;
            window = I(rangeJ, rangeI,:);
            
            % Convert to a vector
            % each column representing a color
            window_vector = reshape(window,win_pixels,3);
            
            % Calculate the mean and difference
            window_mean = mean(window_vector)';
            diff = window_vector' - repmat(window_mean,1,win_pixels);
            
            % both methods of computing covariant produce the same results
            window_covariance = (diff*diff'/win_pixels)+U;
            %window_covariance = (window_vector'*window_vector/win_pixels - window_mean*window_mean')+U;
            
            % Compute the elements in the L matrix
            % easier to just store these in a spares matrix
            L_element = eye(win_pixels) - (1 + diff' * inv(window_covariance) * diff) ./ win_pixels;
            L_element = (L_element(:))'; % reshape it to a single vector
            
            % Calculate the cordinates in the L matrix that we are dealing
            % with.  [ coordinates required in a sparse matrix ]
            
            % Step 1.  Get the indicies of the current window
            window_indicies = reshape(windowIndicies(rangeJ,rangeI),1,win_pixels);
            
            % Step 2.  Create all combinations of pixels
            x = repmat(window_indicies,win_pixels,1);
            y = x';
            
            % reformat combination of pixels
            x = (x(:))';
            y = (y(:))';
            
            
            indicies_x((count*(win_pixels^2) + 1):(count*(win_pixels^2)+(win_pixels^2))) = x;
            indicies_y((count*(win_pixels^2) + 1):(count*(win_pixels^2)+(win_pixels^2))) = y;            
            elements((count*(win_pixels^2) + 1):(count*(win_pixels^2)+(win_pixels^2))) = L_element;
            
            count = count + 1;
        end
    end
    
    L = sparse(indicies_x,indicies_y,elements,dimI(1)*dimI(2),dimI(1)*dimI(2));
    T = (L + .0001 .* speye(size(L))) \ T_est(:) .* .0001;
    T = reshape(T, size(T_est));
end

