%% Haar Wavelet method to solve wave-like partial differential equation

                % KAUSHIK IYER
                
%    y..(x,t) - ((x^2)/2)y''(x,t) = 0   0<x<1 , t>0  
%   where . represents derivative with respect to t an ' represents
%   derivative with respect to x

%  with the initial condition 
%   y(x,0) = f(x) = x,  y.(x,0) = g(x) = x^2
%   y(0,t) = λ_0(t) = 0,    y(1,t) = λ_1(t) = 1 + sinh(t)
             
%           the exact solution is 
%        y(x,t) = x + (x^2)sinh(t) 

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'

% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
x    = (j-0.5) ./ N;                            % grid points
t    = (j(1)-0.5) / N;                                           % initial value
delta_t = 1/N ;                                  % time increments 




% generating the Haar matrix and its integral matrices of size (2M x 2M)
H = zeros(N,N);   % initialising the matrix with zeros
P = zeros(N,N);   % initialising the matrix with zeros

for i = 1:N
    H(i,:) = haar_matrix(x,i,J); 
    P(i,:) = integral_of_H(x,i,J);
    
end
%disp(H);

% defining initial and boundary conditions

E = ones(1,N);

approx_sol = zeros(N,N);

% finding the filter coefficients
for ts = 1:N
    f = zeros(1,N);
    f = (x + t.*E).*cos(x.*t);
    t = t + delta_t;            % time increment 
    
    % finding filter coefficients     
    a(ts,:) = f/(P + delta_t.*H);
    % constructing the approximated solution
    
    approx_sol(ts,:) = delta_t.*(a(ts,:)*P);
                          
end
t = 0;
% true solution
exact_sol = zeros(N,N);
for i = 1:N
    t = t + delta_t;
    for j = 1:N
        exact_sol(i,j) = sin(x(i).*t);
    end    
end

% error
error = zeros(N,N);
error_max  = zeros(1,N);
for i = 1:N
    error(i,:) = abs(exact_sol(i,:) - approx_sol(i,:));
    error_max(i) = max(error(i,:)); % Euclidean norm of the error
end


%% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

figure(1)
surf(approx_sol)
hold on
surf(exact_sol)
xlabel('$x$'); ylabel('$t$'); zlabel('$y(x,t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')


figure(2)
surf(error)
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
            %y(i) = round(2^(j/2)*1,4);
            y(i) = 1;
        elseif (beta <= t(i) && (t(i) < gamma))
            %y(i) = round(2^(j/2)*(-1),4);
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







function F = myfun(a,x,Q,Q_1,lambda_0,lambda_1_d2_tn1,y_d2,E)
F = a*(Q - x.*Q_1) - (((x.^2).*y_d2)./2 + x.*(lambda_0 - lambda_1_d2_tn1) - lambda_0.*E);

end





