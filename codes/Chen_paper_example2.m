%% Haar Wavelet method to solve fractional differential equations

                % KAUSHIK IYER
                
%  D^(1/4)u(t) + u(t) = 1.5*t^(2/3)/Γ(2/3) + t^(4/3)           
% with the initial condition u(0) = 0          
             
%           the exact solution is 
%        u(t) = t

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (j-0.5) ./ N;                                % grid points
alpha = input('Enter the value of alpha:');


%generating the Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
for i = 1:N
    H(i,:) = haar_matrix(t,i,J); 
end
%disp(H);

% generating the F_α matrix (NxN)
F_alpha_initial = diag(ones(1,N));
zeta = zeros(1,N);
for i = 1:N
    zeta(i) = (i+1)^(alpha+1) - 2*(i^(alpha+1)) + (i-1)^(alpha+1);
end

for row = 1:N
    for col = row+1:N
        F_alpha_initial(row,col) = zeta(col-row); 
    end
end
F_alpha = (1/((N^alpha)*gamma(alpha+2))).*F_alpha_initial;

% generating the operational matrix (NxN)
P = zeros(N,N);
P_alpha = H*F_alpha/(H);



% finding the filter coefficients
A = diag(ones(1,N));                             % coefficient function matrix at the collocation points
f = t.^4 - 0.5*t.^(3) - (3/(gamma(4-alpha))).*t.^(3-alpha) + (24/(gamma(5-alpha))).*t.^(4-alpha);    % function values at the collocation points
c = f/(H + P_alpha*H*A);

% constructing the approximated solution 
approx_sol = c*P_alpha*H;

% true solution
exact_sol = t.^4 - 0.5*t.^3;

% error
error = abs(exact_sol - approx_sol);
error_norm = norm(error,2); % Euclidean norm of the error



%% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

oft = 0.01;

% fig:01
figure('color','w')
plot(t,approx_sol,'g',t,exact_sol,'rs')
xlabel('$t$'); ylabel('$u(t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')
axis([-oft 1+oft min(exact_sol)-oft max(exact_sol)+oft])


% fig:02
figure('color','w')
plot(t,error,'r.-')
xlabel('$t$'); ylabel('Absolute Error');
title('Absolute Error: $\max|u_{numeric} - u_{analytic}|$')
axis([-oft 1+oft min([error])-oft max([error])+oft])




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
            y(i) = round(2^(j/2)*1,4);
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = round(2^(j/2)*(-1),4);
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

function F = myfun(a,H,Q,P,E)
F = [a*(H+Q)+E];
end

function y = q_1(t,i,J,A,B)
% Function to generate the Q_1 function for the given interval for a
% given index 'i'
% It is defined as            [0.5(B-A)^2 , i=1  
%                     Q_i(1) =[
%                             [ 0.25*(B-A)^2/(m^2) , i>1

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the Q_1 function

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
       y(i) = 0.5*(B-A)^2;     
    end
else
    for i = 1:length(t)
         y(i) = (B-A)^2 / (4*m^2);
    end
end
end




