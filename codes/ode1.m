%% Haar Wavelet method to solve fractional partial differential equations

                % KAUSHIK IYER
                
%  u'(t) - u(t) = exp(t)         
% with the initial condition u(0) = 0          
             
%           the exact solution is 
%        u(t) = t*exp(t)
clear;
clc;
addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (j-0.5)./N;
y0 = 0;
%generating the Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
P = H;
for i = 1:N
    P(i,:) = integral_of_H(t,i,J);
    H(i,:) = haar_matrix(t,i,J);
end


approx_sol = zeros(1,length(t));
exact_sol  = zeros(1,length(t)); 
error      = zeros(1,length(t));
E          = zeros(1,length(t));
c1         = zeros(1,length(t));

        
% true solution
for i = 1:N
    exact_sol(i) = t(i)*exp(t(i));
end

% finding filter coefficients
c1 = (exp(t))/(H - P);

% constructing the approximated solution 
approx_sol = c1*P;    
   

% error
error = abs(exact_sol - approx_sol);    
   


% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

% fig:01
figure('color','w')
plot(t,approx_sol,'rs',t,exact_sol,'g')
xlabel('$t$');ylabel('$y(t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')



% fig:02
figure('color','w')
plot(error)
xlabel('$t$');ylabel('Absolute Error');
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


