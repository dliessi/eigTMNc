function [dX, dY, tau, AXY, AYY, BXY, BYY, CXX, CXY, CYX, CYY] = ...
         modelc_expercoup(~)
% modelc_expercoup defines a coupled RE/RFDE model for eigTMNc.
%
% It defines the model
% x(t) = - 1/2 [int_{0}^{7/2 pi} x(t - theta) d theta
%              - int_{0}^{pi/2} ln(y(t - theta)) d theta]
% y'(t) = - ln(y(t - pi/2)) y(t)
% linearized around its periodic solution (sin(t), exp(sin(t)) .

% Lines not marked with '% INPUT' should not be changed.


% Mnemonics for parameters.
% par(1) = ...
% par(2) = ...
% ...

% eigTMNc passes the par variable to modelc.
% In case it is necessary because the model definition (not just the
% values of the A**, B** and C** matrices) depends on the parameter
% values, you can replace the ~ in the first line with a variable
% name and use it inside the modelc function.


% Model dimensions: RE...
dX = 1; % INPUT
% ... and RFDE.
dY = 1; % INPUT

% Delays (includes 0 and is sorted: [0 tau_{1} tau_{2} ...]).
tau = [0 pi/2 7*pi/2]; % INPUT


% Number of nontrivial delays.
p = length(tau) - 1;

% Initialize A** matrices as empty matrices.
AXY = [];
AYY = [];

% Initialize cell arrays for B** and C** matrices as empty matrices.
BXY = cell(p, 1);
BYY = cell(p, 1);
CXX = cell(p, 1);
CXY = cell(p, 1);
CYX = cell(p, 1);
CYY = cell(p, 1);

% Define the model.
% A**'s and B**{*}'s must be function handles taking (t, par) as
% arguments.
% C**{*}'s must be function handles taking (t, theta, par) as
% arguments.

% AXY = @(t, par) []; % INPUT

AYY = @AYY; % INPUT

% BXY{1} = @(t, par) []; % INPUT
% ...
% BXY{p} = @(t, par) []; % INPUT

BYY{1} = @BYY1; % INPUT
% ...
% BYY{p} = @(t, par) []; % INPUT

CXX{1} = @CXX1; % INPUT
CXX{2} = @CXX2; % INPUT
% ...
% CXX{p} = @(t, theta, par) []; % INPUT

CXY{1} = @CXY1; % INPUT
% ...
% CXY{p} = @(t, theta, par) []; % INPUT

% CYX{1} = @(t, theta, par) []; % INPUT
% ...
% CYX{p} = @(t, theta, par) []; % INPUT

% CYY{1} = @(t, theta, par) []; % INPUT
% ...
% CYY{p} = @(t, theta, par) []; % INPUT

% For performance reasons you may want to define subfunctions and
% pass their function handle instead of using anonymous functions.

end


function AYY = AYY(t, ~)
AYY = - log(exp(sin(t - pi/2)));
end

function BYY1 = BYY1(t, ~)
BYY1 = - exp(sin(t)) / exp(sin(t - pi/2));
end

function CXX1 = CXX1(~, ~, ~)
CXX1 = - 0.5;
end

function CXX2 = CXX2(~, ~, ~)
CXX2 = - 0.5;
end

function CXY1 = CXY1(t, theta, ~)
CXY1 = 0.5 / exp(sin(t + theta));
end
