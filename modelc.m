function [dX, dY, tau, AXY, AYY, BXY, BYY, CXX, CXY, CYX, CYY] = ...
         modelc(~)
% modelc defines a coupled RE/RFDE model for eigTMNc.
%
% It defines the model
% x(t) = sum_{k = 1}^{p}
%            int_{- tau_{k}}^{- tau_{k - 1}}
%                C_{XX}^{(k)}(t, theta) * x(t + theta) d theta
%        + A_{XY}(t) * y(t)
%        + sum_{k = 1}^{p} B_{XY}^{(k)}(t) * y(t - tau_{k})
%        + sum_{k = 1}^{p}
%              int_{- tau_{k}}^{- tau_{k - 1}}
%                  C_{XY}^{(k)}(t, theta) * y(t + theta) d theta
% y(t) = sum_{k = 1}^{p}
%            int_{- tau_{k}}^{- tau_{k - 1}}
%                C_{YX}^{(k)}(t, theta) * x(t + theta) d theta
%        + A_{YY}(t) * y(t)
%        + sum_{k = 1}^{p} B_{YY}^{(k)}(t) * y(t - tau_{k})
%        + sum_{k = 1}^{p}
%              int_{- tau_{k}}^{- tau_{k - 1}}
%                  C_{YY}^{(k)}(t, theta) * y(t + theta) d theta

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
tau = [0 1 3]; % INPUT


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

% AYY = @(t, par) []; % INPUT

% BXY{1} = @(t, par) []; % INPUT
% ...
% BXY{p} = @(t, par) []; % INPUT

% BYY{1} = @(t, par) []; % INPUT
% ...
% BYY{p} = @(t, par) []; % INPUT

% CXX{1} = @(t, theta, par) []; % INPUT
% ...
% CXX{p} = @(t, theta, par) []; % INPUT

% CXY{1} = @(t, theta, par) []; % INPUT
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
