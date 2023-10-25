%% Haar Wavelet method to solve integro differenial equations

                % KAUSHIK IYER
                
%           u(x) - integral[K(x,t)*u(t)dt] = f(x)           
% here K(x,t) = x + t
% f(x) = x^2;
             
%           the exact solution is 
%        u(x) = x^2 - 5x-17/6

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J+1); % N = 2M                              % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (j-0.5) ./ N;                                % grid points


%generating the Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
for i = 1:N
    H(i,:) = haar_matrix(t,i,J); 
end
%disp(H);

% generating the G matrix (2M x 2M)
G = zeros(N,N);
for i = 1:N
    G(i,:) = integral_of_H_at_1(t,i,J);
end


% finding the filter coefficients
a = zeros(1,length(t));             % coefficient matrix definition
f = t.^2;                           % function values at the collocation points

a = f/(H-G);                      

% constructing the approximated solution 
approx_sol = a*H;

% true solution
exact_sol = t.^2 - 5.*t - 17/6;

% error
error = zeros(1,length(t));
for i = 1:length(t)
    error(1,i) = abs(exact_sol(i) - approx_sol(i));
    error_norm = norm(error,2); % Euclidean norm of the error
end 
max_error = max(error);

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
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')
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
%           y(i) = round(2^(j/2)*1,4);
            y(i) = 1;
        elseif (beta <= t(i) && (t(i) < gamma))
%           y(i) = round(2^(j/2)*(-1),4);
            y(i) = -1;
        else
            y(i) = 0;
        end
    end
end


end


function y = integral_of_H_at_1(t,i,J)
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

    [ind_j, ind_i] = find(IMask == i);
    m = 2^(ind_j - 1);
    k = ind_i - 1;
end


y = zeros([1 length(t)]);

if i == 1  
    for j = 1:length(t)
        y(i,j) = t(j) + 0.5;
    end
 else  
     for i = 1:length(t)
         y(i) = -1/(4*(m^2));
     end
end
end


