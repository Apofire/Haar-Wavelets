%% Haar Wavelet method to solve ODE

                % KAUSHIK IYER
                
% A second order ODE hasbeen considered here 
%            y'' + y' = 0    
% with the initial conditions y(0) = 1 ,y'(0) = 0
 
%           the exact solution is 
%           y(x) = cos(x)

%         ------method of solving------
%         y'' = aH
%         y'  = aP + y'(0)
%         y   = aQ + y'(0)x + y(0)  
% then substitute these equations in the ODE and solve for the filter coefficients
addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J = 3;                              % level of decomposition
N = 2^(J + 1); % N = 2M             % number of basis functions
j = 1:N;                            % index of grid points
x = (j - 0.5) ./ N;                 % grid points

% initialising the initial conditions of the ODE
y_0 = 1;                            % y(0) = 1
y1_0 = 0;                           % y'(0) = 0

%generating the Haar matrix of size (2M x 2M)
H = zeros(N,N);                     % initialising the matrix with zeros

for i = 1:N
    H(i,:) = haar_matrix(x,i,J); 
end
%disp(H);

%generating the integral of the Haar matrix (2M x 2M)
P = zeros(N,N);

for i = 1:N
    P(i,:) = integral_of_H(x,i,J); 
end
%disp(P);

%generating the integral of the Haar matrix (2M x 2M)
Q = zeros(N,N);

for i = 1:N
    Q(i,:) = integral_of_P(x,i,J); 
end
%disp(Q);

% finding the filter coefficients
E    = ones(1,N);
%a    = fsolve(@(a) myfun(x,H,P,Q,E),0.5.*ones(N,1));
a    = -E/(H+Q);

y1 = sin(x);
approx = (sin(x)/H)*Q;
exact = -sin(x) + x.*ones(1,N);

% constructing the approximated solution
y = zeros(N,1);
for j = 1:N    
    S = 0;
    for i = 1:N
        S = S + a(i) * integral_of_P(x(j),i,J);
    end
    y(j) = y_0 + x(j) * y1_0 + S;
end % for

% exact solution
y_exact = cos(x); 

% %% Runge - Kutta method
% [x, y1] = ode113('model1', x, [y_0 y1_0]);

%% Plot graphics
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

oft = 0.01;

% fig:01
figure('color','w')
plot(x,y_exact,'g',x,y,'rs')
xlabel('$t$'); ylabel('$y$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Exact','WHM')
axis([-oft 1+oft min(y_exact)-oft max(y_exact)+oft])

% Absolute errors
rWHM = abs(y - y_exact');
% fig:02
figure('color','w')
plot(x,rWHM,'r.-')
xlabel('$t$'); ylabel('Absolute Error');
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')
legend('WHM',...
    'Location','northoutside','Orientation','horizontal')


%% Disp Errors
disp(['error RGK: ' num2str(max(rRGK)) ' error WHM: ' num2str(max(rWHM)) ...
    ' error RW: ' num2str(max(rRW))])

%% Functions

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
            y(i) = 1;
        elseif (beta <= x(i) && (x(i) < gamma))
            y(i) = -1;
        else
            y(i) = 0;
        end
    end
end
%disp(y);
end


function y = integral_of_H(x,i,J)
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


y = zeros([1 length(x)]);

if i == 1    
    for i = 1:length(x)
        if (0 <= x(i) && (x(i) < 1))
            y(i) = x(i);
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
            y(i) = x(i) - alpha;
        elseif (beta <= x(i) && (x(i) < gamma))
            y(i) = gamma - x(i);
        else
            y(i) = 0;
        end
    end
end
end


function y = integral_of_P(x,i,J)
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


y = zeros([1 length(x)]);

if i == 1    
    for i = 1:length(x)
        if (0 <= x(i) && (x(i) < 1))
            y(i) = 0.5 * x(i) * x(i);
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
            y(i) = 0.5 * (x(i) - alpha).^2;
        elseif (beta <= x(i) && (x(i) < gamma))
            y(i) = 1 / (4*m^2) - 0.5 * (gamma - x(i)).^2;
        elseif (gamma <= x(i) && (x(i) < 1))
            y(i) = 1 / (4*m^2);
        else
            y(i) = 0;
        end
    end
end
end

function F = myfun(a,H,Q,P,E)
F = [a*(H+Q)+E];
end

