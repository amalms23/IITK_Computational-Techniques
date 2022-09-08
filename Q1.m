%##############################################
%#############  NEWTON RAPHSON   ##############
%##############################################

f = @(x) x^3-0.5;     % function definition
derf = @(x) 3*x^2;    % derivative of the function

tic
x0 = 0;              % initial guess 
i = 1;               % counter
er = 10;             % error for stop criterion


while er > 10^-6 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #
    
    i = i+1;                              % counter
    x1 = x0 - f(x0)/derf(x0);       % N-R Algorithm
    er = abs(x1 - x0);          % error calculation
    x0 = x1;          % updating the initial values
end


if (er > 10^-6) || isnan(x1)
    fprintf("For intial approximation = 0, " + ...
        "Newton Raphson method for finding" + ...
        " the root **FAILED** \n");
else
    fprintf("For intial approximation = 0, " + ...
        "root using NR Method is = %f and no: of iterations =" + ...
        " %d\n",x1,i);
end
time = toc


tic

x0 = 1;              % initial guess 
i = 1;               % counter
er = 10;             % error for stop criterion


while er > 10^-6 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #
    
    i = i+1;                              % counter
    x1 = x0 - f(x0)/derf(x0);       % N-R Algorithm
    er = abs(x1 - x0);          % error calculation
    x0 = x1;          % updating the initial values
end


if (er > 10^-6) || isnan(x1)
    fprintf("For intial approximation = 1, " + ...
        "Newton Raphson method for finding" + ...
        " the root **FAILED** \n");
else
    fprintf("For intial approximation = 1, " + ...
        "root using NR Method is = %f and no: of iterations =" + ...
        " %d\n",x1,i);
end
time = toc

tic

x0 = 0.5;              % initial guess 
i = 1;               % counter
er = 10;             % error for stop criterion


while er > 10^-6 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #
    
    i = i+1;                              % counter
    x1 = x0 - f(x0)/derf(x0);       % N-R Algorithm
    er = abs(x1 - x0);          % error calculation
    x0 = x1;          % updating the initial values
end


if (er > 10^-6) || isnan(x1)
    fprintf("For intial approximation = 0.5, " + ...
        "Newton Raphson method for finding" + ...
        " the root **FAILED** \n");
else
    fprintf("For intial approximation = 0.5, " + ...
        "root using NR Method is = %f and no: of iterations =" + ...
        " %d\n",x1,i);
end

time = toc

%##############################################
%#############   MULLER METHOD   ##############
%##############################################


f = @(x) x^3-0.5;     % function definition

xk20 = 0;            % initial guess 
xk10 = 1;            % initial guess
xk0 = 0.5;

xk2 = 0;            % initial guess 
xk1 = 1;            % initial guess
xk = 0.5;          % initial guess
i = 1;               % counter
er = 10;             % error for stop criterion

tic
while er > 10^-6 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #

    i = i+1;                              % counter
    lk = (xk-xk1)/(xk1-xk2);
    delk = lk + 1;
    gk = (((lk)^2)*f(xk2))-((delk)^2)*f(xk1)+(lk+delk)*f(xk);
    ck = lk*((lk*f(xk2))-((delk)*f(xk1))+f(xk));
    if gk>0
    lk1 = (-2*delk*f(xk))/(gk+sqrt((gk^2)-4*delk*ck*f(xk)));
    else
    lk1 = (-2*delk*f(xk))/(gk-sqrt((gk^2)-4*delk*ck*f(xk)));
    end
    xk11 = xk+(xk-xk1)*lk1;

    er = abs(xk11 - xk);          % error calculation

        xk2 = xk1;
        xk1 = xk;
        xk = xk11;          % updating the initial values
end


if (er > 10^-6) || isnan(xk11)

    fprintf("For intial approximation = %d, %d, & %d " + ...
        "Mullers method for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xk20,xk10,xk0,i);

else

    fprintf("For the given intial approximation = %.2f, %.2f, & %.2f, " + ...
        "root = %f and no: of iterations =" + ...
        " %d\n",xk20,xk10,xk0,xk11,i);

end
time = toc


%##############################################
%#############    Sixth Order    ##############
%##############################################


f = @(x) x^3-0.5;     % function definition
fd = @(x) 3*x^2;

xn0 = 0;
xn = 0;          % initial guess
i = 1;               % counter
er = 10;             % error for stop criterion

tic
while er > 10^-9 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #

    i = i+1;                              % counter
    
    wn = xn - f(xn)/fd(xn);
    zn = wn-(f(wn)/fd(xn))*((2*f(xn)-f(wn))/(2*f(xn)-5*f(wn)));
    xn1 = zn-(f(zn)/fd(xn))*((f(xn)-f(wn))/(f(xn)-3*f(wn)));

    er = abs(xn1 - xn);          % error calculation

        
        xn = xn1;          % updating the initial values
end


if (er > 10^-4) || isnan(xn1)

    fprintf("For intial approximation = %d, " + ...
        "sixth order method for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xn0,i);

else

    fprintf("For intial approximation = %d, " + ...
        "root using sixth order method is = %f and no: of iterations =" + ...
        " %d\n",xn0,xn1,i);

end
time = toc

xn0 = 1;
xn = 1;          % initial guess
i = 1;               % counter
er = 10;             % error for stop criterion

tic
while er > 10^-9 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #

    i = i+1;                              % counter
    
    wn = xn - f(xn)/fd(xn);
    zn = wn-(f(wn)/fd(xn))*((2*f(xn)-f(wn))/(2*f(xn)-5*f(wn)));
    xn1 = zn-(f(zn)/fd(xn))*((f(xn)-f(wn))/(f(xn)-3*f(wn)));

    er = abs(xn1 - xn);          % error calculation

    xn = xn1;          % updating the initial values
end


if (er > 10^-4) || isnan(xn1)

    fprintf("For intial approximation = %d, " + ...
        "sixth order for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xn0,i);

else

    fprintf("For intial approximation = %d, " + ...
        "root using sixth order method is = %f and no: of iterations =" + ...
        " %d\n",xn0,xn1,i);

end
time = toc

xn0 = 1/2;
xn = 1/2;          % initial guess
i = 1;               % counter
er = 10;             % error for stop criterion
tic

while er > 10^-9 & i < 100 
%(the condition for stopping the loop is either: #
%   1. the error(er) < 10^-6                 or  #
%   2. the total number of iterations >=100      #

    i = i+1;                              % counter
    
    wn = xn - f(xn)/fd(xn);
    zn = wn-(f(wn)/fd(xn))*((2*f(xn)-f(wn))/(2*f(xn)-5*f(wn)));
    xn1 = zn-(f(zn)/fd(xn))*((f(xn)-f(wn))/(f(xn)-3*f(wn)));

    er = abs(xn1 - xn);          % error calculation
    xn = xn1;          % updating the initial values
end


if (er > 10^-4) || isnan(xn1)

    fprintf("For intial approximation = %d, " + ...
        "sixth order method for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xn0,i);

else

    fprintf("For intial approximation = %.2f, " + ...
        "root using sixth order method is = %f and no: of iterations =" + ...
        " %d\n",xn0,xn1,i);

end
time = toc

