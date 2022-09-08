%#################################################%
%############# NEWTON-RAPHSON METHOD #############%
%#################################################%

f = @(x) log10(x)-x+3;     % function definition
derf = @(x) (1/((log(10))*x))-1;    % derivative of the function

%#################################################%

x0 = 0.25;              % initial guess 
i = 1;               % counter
er = 10;             % error for stop criterion

tic
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
    fprintf("For intial approximation = 0.25, " + ...
        "Newton Raphson method for finding" + ...
        " the root **FAILED** \n");
else
    fprintf("For intial approximation = 0.25," + ...
        "root using Newton Raphson method is = %f and no: of iterations =" + ...
        " %d\n",x1,i);
end
time = toc
%#################################################%

x0 = 0.5;              % initial guess 
i = 1;               % counter
er = 10;             % error for stop criterion

tic
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
        "root using Newton Raphson method is = %f and no: of iterations =" + ...
        " %d\n",x1,i);
end
time = toc
%#################################################%

x0 = 1;              % initial guess 
i = 1;               % counter
er = 10;             % error for stop criterion

tic
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
        "root using Newton Raphson method is = %f and no: of iterations =" + ...
        " %d\n",x1,i);
end
time = toc
%#################################################%
%#############     Muller METHOD     #############%
%#################################################%

f = @(x) log10(x)-x+3;     % function definition

xk20 = 0.25;            % initial guess 
xk10 = 0.5;            % initial guess
xk0 = 1;

xk2 = 0.25;            % initial guess 
xk1 = 0.5;            % initial guess
xk = 1;          % initial guess
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
        "root using Muller method is = %f and no: of iterations =" + ...
        " %d\n",xk20,xk10,xk0,xk11,i);

end
time = toc
%#################################################%
%#############   Sixth Order METHOD  #############%
%#################################################%


f = @(x) log10(x)-x+3;     % function definition
fd = @(x) (1/(log(10)*x))-1;

xn0 = 0.25;
xn = 0.25;          % initial guess
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

    fprintf("For intial approximation = %.2f, " + ...
        "sixth order method for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xn0,i);

else

    fprintf("For intial approximation = %.2f, " + ...
        "root using sixth order method = %f and no: of iterations =" + ...
        " %d\n",xn0,xn1,i);

end
time = toc

xn0 = 0.5;
xn = 0.5;          % initial guess
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

    fprintf("For intial approximation = %.2f, " + ...
        "sixth order for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xn0,i);

else

    fprintf("For intial approximation = %.2f, " + ...
        "root using sixth order method = %f and no: of iterations =" + ...
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
        "sixth order method for finding" + ...
        " the root **FAILED** at iteration number: %d\n",xn0,i);

else

    fprintf("For intial approximation = %.2f, " + ...
        "root using sixth order method is = %f and no: of iterations =" + ...
        " %d\n",xn0,xn1,i);

end
time = toc

