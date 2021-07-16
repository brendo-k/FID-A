% Rosette.m
% Brenden Kadota, McGill University 2019.
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
% par.sw = spectral width
% par.Fov = Fov array in m, should be 2x1 dimension
% par.ImgeSize = Number of points in the x and y dimensions, 2x1 dimension.
%
% OUTPUTS:
% obj = Trajectory obj

function [obj] = Rosette(par, file_path)

    %some example numbers taken from Schirda et al. 2019
    if(nargin < 1)
        par.sw = 2500;
        par.Fov = [0.2,0.2];
        par.imageSize = [16, 16, 1024];
    end
    
    %initalize constant
    gamma = 42.577478518e6;    
    
    %get omega1 and omega2 from sw
    omega1 = par.sw*pi;
    omega2 = omega1;
    %get kFov
    kFov = par.imageSize(1)/par.Fov(1);
    %get KMax
    kMax = kFov/2;
    %calculate gMax
    gMax = kMax * max(omega1, omega2) / (gamma); %T/m
  
    %dwell time calculation
    dwellTime = 1/(gamma*gMax*par.Fov(1) + par.sw);
    %number of angular interleves (number of turns)
    N_AngleInts = floor((pi*par.imageSize(1))/(sqrt(1+3*(omega2^2/omega1^2))));
    
    %get time vector
    t = 0:dwellTime:(1/par.sw)*(par.imageSize(3)-1);
    
    %initalize array size
    traj = complex(zeros(N_AngleInts, length(t)), 0);
    %Rosette trajectory formulat from shirda paper
    traj(1,:) = kMax*sin(omega1*t).*exp(1i*omega2*t);

    %creating rosette trajectory in k space
    for rotation = 1:N_AngleInts
        %get the radians to rotate each trajectory by
        rotationAngle = ((rotation - 1)/N_AngleInts)*2*pi;
        
        %rotate trajecotry in the complex plane using e^(i*theta)
        rotationConverter = exp(1i*rotationAngle); 
        
        %apply to the first trajectory to get the rotated trajectory.
        traj(rotation,:) = rotationConverter*traj(1,:);
    end
    if(exist('file_path', 'var'))
        svg = zeros(numel(traj),3);
        writematrix(['TR', 'K Space', 'time'], file_path, 'WriteMode', 'overwrite', 'Delimiter',',');
        
        linear_traj = traj';
        linear_traj = linear_traj(:);
        tr = repmat([1:size(traj,1)], [size(traj,2), 1]);
        tr = tr(:);
        time = repmat([0:dwellTime:dwellTime*(size(traj,2)-1)]', [size(traj,1),1]);
        svg = horzcat(tr, linear_traj, time);
        writematrix(svg, file_path, 'WriteMode', 'overwrite', 'Delimiter', ',');
        
        
       
        writematrix(svg, file_path, 'WriteMode', 'append');
                
    
    end
    
    %Create trajectory object
    obj = Trajectory(traj, par.imageSize, par.Fov, dwellTime, ['x', 'y'], par.sw);
end
