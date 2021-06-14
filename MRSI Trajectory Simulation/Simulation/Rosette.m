% Rosette.m
% Jamie Near, McGill University 2019.
% 
% USAGE:
% Rosette(par)
% 
% DESCRIPTION:
% This function takes in a parameter sturucture (par) and outputs the
% corresponding rosette trajectory form the parameters. The formulas to
% create the rosette trajectory comes from Schirda et al. 2009.
%
% INPUTS:
% par = parameter structure
% par.omega1 = angular turning rate of trajectory
% par.omega2 = osilitory rate of trajectroy. (how fast the trajectory
% completes one cicle in the rosette)
% par.kMax = maxium k space coordinates
% par.dwellTime = dwell time of scan
% par.readOutTime = read out time
% par.cycleTime = one spectral dwell time. (ie. the time to complete one
% cycle)
% par.slewRate = gradient slew rate
% par.nAngInts = number of angular interleaves
% par.repetitionTime = TR, time between scans.
%
% OUTPUTS:
% inputTraj   = the original k space trajectory
% gradentTraj     = the gradient trajectory
% finalKSpaceTraj = the trajectory after acounting for scanner slew rates
% gMax = the max gradient needed for the trajectory
% maxSlewRate = the max slew rate needed for the trajectory (excluding the
% inital ramping).

function [obj] = Rosette(par)

    %some example numbers taken from Schirda et al. 2019
    if(nargin < 1)
        par.sw = 2500;
        par.Fov = [0.2,0.2];
        par.imageSize = [16, 16, 1024];
    end
    
    %initalize constant
    gamma = 42.577478518e6;    
    
    omega1 = par.sw*pi;
    omega2 = omega1;
    kFov = par.imageSize(1)/par.Fov(1);
    kMax = kFov/2;
    gMax = kMax * max(omega1, omega2) / (gamma); %T/m
  
    dwellTime = 1/(gamma*gMax*par.Fov(1) + par.sw);
    N_AngleInts = floor((pi*par.imageSize(1))/(sqrt(1+3*(omega2^2/omega1^2))));
    
    t = 0:dwellTime:(1/par.sw)*(par.imageSize(3)-1);
    
    traj = complex(zeros(N_AngleInts, length(t)), 0);
    traj(1,:) = kMax*sin(omega1*t).*exp(1i*omega2*t);

    %creating rosette trajectory in k space
    for rotation = 1:N_AngleInts
        rotationAngle = ((rotation - 1)/N_AngleInts)*2*pi;
        rotationConverter = exp(1i*rotationAngle); 
        traj(rotation,:) = rotationConverter*traj(1,:);
    end
    
    obj = Trajectory(traj, par.imageSize, par.Fov, dwellTime, ['x', 'y'], par.sw);
end
