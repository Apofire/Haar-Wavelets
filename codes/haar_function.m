J    = 1;                                  % level of decomposition
N    = 2^(J+1); % N = 2M                   % number of basis functions
j    = 1:N;                                % index of grid points
x    = (j-0.5) ./ N;                       % grid points

H = zeros(N,N);                     % initialising the matrix with zeros
for i = 1:N
    H(i,:) = haar_matrix(x,i,J); 
end
disp(H);

function y = haar_matrix(x,i,J)
% Function to generate the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the haar function

if i == 1
    m = 0;
    k = 0;
    j = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [ind_j, ind_i] = find(IMask == i);
    m = 2^(ind_j - 1);
    k = ind_i - 1;
    j = ind_j-1;
end

y = zeros([length(x) 1]);

if i == 1    
    for i = 1:length(x)
        if (0 <= x(i) && (x(i) < 1))
            y(i) = 1;
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(x)
        if (alpha <= x(i) && (x(i) < beta))
            y(i) = round((2^(j/2))*1,4);
        elseif (beta <= x(i) && (x(i) < gamma))
            y(i) = round((2^(j/2))*(-1),4);
        else
            y(i) = 0;
        end
    end
end


end