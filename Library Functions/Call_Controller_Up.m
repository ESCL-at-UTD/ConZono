function [ Output ] = Call_Controller_Up( Input, plot, time )

Output = Input;

%% Construct Decision Variable Input to Controller

Inputs(1) = {Output.x0};
Inputs(2:Output.N+1) = mat2cell(Output.rdes,Output.Nr,ones(1,Output.N));
Inputs(Output.N+2:2*Output.N+2) = mat2cell(Output.t,1,ones(1,Output.N+1));
Inputs(2*Output.N+3:3*Output.N+3) = mat2cell(Output.tf,1,ones(1,Output.N+1));

%% Calculate Decision Variable Outputs from Controller
% Solve the optimization problem
if time == 1
    tic;
    [OUT,diagnostics] = Output.Controller{Inputs};
    Output.t_calc = toc;
    if diagnostics == 0
        display([Output.Name, ' controller calc time is ', num2str(Output.t_calc), ' seconds.'])
    else
        display([Output.Name, ' Yalmip error ', num2str(diagnostics)])
        Output.t_calc = 0;
    end
else
    OUT = Output.Controller{Inputs};
end

Output.x = cell2mat(OUT(1:Output.N+1));
Output.u = cell2mat(OUT(Output.N+2:2*Output.N+1));
                        
%% Plot Results

if plot == 1
    figure;
    subplot(3,2,1);stairs(0:Output.N,Output.x(1,1:end));
    subplot(3,2,3);stairs(0:Output.N,Output.x(2,1:end));
    subplot(3,2,5);stairs(0:Output.N,Output.x(3,1:end));
    subplot(3,2,2);stairs(Output.u(1,:));
    subplot(3,2,4);stairs(Output.u(2,:));
    subplot(3,2,6);stairs(Output.rdes(end,:)-Output.u(3,:));

end
