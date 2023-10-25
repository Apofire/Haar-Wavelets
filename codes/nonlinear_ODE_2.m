%% Haar Wavelet method to solve system of nonlinear ordinary differential equations

                % KAUSHIK IYER
                
%           y1''(t) + 2/t*y'(t) + (y(t)^n)= 0, y(0) = 1
%           y'(0) = 0,  
             
%           the exact solution is 
%        y(t) =  exp(-t^2)    

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J+1); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (2*j-1) ./ (2*N);                                % grid points



%generating the Haar matrix and integral of Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
P = zeros(N,N);
Q = zeros(N,N);
for i = 1:N
    H(i,:) = haar_matrix(t,i,J); 
    P(i,:) = integral_of_H(t,i,J);
    Q(i,:) = integral_of_P(t,i,J);
end
%disp(H);



% finding the filter coefficients
E  = ones(1,length(t)); 

a = fsolve(@(a) myfun(a,t,H,P,Q,E), zeros(1,N));


% constructing the approximated solution 
approx_sol = a*Q + E;


% true solution
exact_sol = (E + (t.^2)./3).^(-1/2);

% error
error = zeros(length(t),1);
error = abs(exact_sol - approx_sol);


% Euclidean norm of the error
error_norm = norm(error,2);
 




 %% Plot graphics 
% fig:01
figure('color','w')
plot(t,approx_sol,'g',t,exact_sol,'rs')
xlabel('$t$'); ylabel('$y(t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')


% fig:02
figure('color','w')
plot(t,error,'r.-')
xlabel('$t$'); ylabel('Absolute Error');
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')





%% Functions

function y = haar_matrix(t,i,J)
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
    j = ind_j -1;
end

y = zeros([length(t) 1]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = 1;
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = 1;
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = -1;
        else
            y(i) = 0;
        end
    end
end


end


function y = integral_of_H(t,i,J)
% Function to generate the integral of the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the integral of the haar function

if i == 1
    m = 0;
    k = 0;
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

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end


y = zeros([1 length(t)]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = t(i);
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = t(i) - alpha;
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = gamma - t(i);
        else
            y(i) = 0;
        end
    end
end



end

function y = integral_of_P(t,i,J)
% Function to generate the integral of the P function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the integral of the P function

if i == 1
    m = 0;
    k = 0;
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

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end


y = zeros([1 length(t)]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = 0.5 * t(i) * t(i);
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = vpa(0.5 * (t(i) - alpha).^2);
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = vpa(1 / (4*m^2) - 0.5 * (gamma - t(i)).^2);
        elseif (gamma <= t(i) && (t(i) < 1))
            y(i) = vpa(1 / (4*m^2));
        else
            y(i) = 0;
        end
    end
end
end

function F = myfun(a,t,H,P,Q,E)
F = a*H + (2./t).*(a*P) + ((E + a*Q).^5);
end


