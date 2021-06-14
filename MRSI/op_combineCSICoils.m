%op_combineCSICoils.m
%
% Combines the coils of MRSI data. First, coil's phases are adjusted to
% match the phase of the first coil. Then, sensetivity weights for coils are
% calculated by the formula w(i) = S(i)/sqrt(sum(S^2)) where S is intensity of the
% first point of the fids in the center of k space and i is the coil
% number. Coils are then summed by F = sum(W(i)*f(i)) where W(i) is coil's
% sensetivity weight at coil i, f(i) is the fids at coil i and F is the
% final summed signal.
% 
% USAGE:
% in = op_combineCSICoils(in, (optional) samplePoint)
% 
% INPUT:    
% in            = input MRSI object for coils to be combined
% samplePoint   = samplePoint from fid to be used for phase
%               correction. Default to first point.
% samplePoint   = point in in.fids to be used for phase correction
%               of the coils 
%
% OUTPUT:
% in            = MRSI object with combined coils in fids and specs.

function [in, phases, weights] = op_combineCSICoils(in, phaseMap, weightMap)

            
    %some pre condition checks and setting default values
    
    if in.flags.addedrcvrs == 1
        error('coils already combined!')
    end
    
    if in.flags.spatialFT ~= 1 || in.flags.spectralFT ~= 0
        error('Please us op_CSIFourierTransform along the spatial dimension')
    end
    
    samplePoint = 1;
    
    if ~exist('phaseMap', 'var')
        %phase map arranged by coils, x coordinate, y coordinate
        phaseArray = zeros([in.sz(in.dims.coils), in.sz(in.dims.x), in.sz(in.dims.y)]);
        
        for x = 1:in.sz(in.dims.x)
            for y = 1:in.sz(in.dims.y)
                for i = 1:in.sz(in.dims.coils)
                    %calculate the phase of each each coordinate for all the
                    %coils
                    phaseArray(i,x,y) = phase(in.specs(samplePoint,i,x,y));
                end
            end
        end
        
        phases = phaseArray; 

        phaseArray = repmat(phaseArray, [1, 1, 1, in.sz(in.dims.t)]);
        phaseArray = permute(phaseArray, [4, 1, 2, 3]);
    else       
        phaseArray = repmat(phaseMap, [1, 1, 1, in.sz(in.dims.t)]);
        phaseArray = permute(phaseArray, [4, 1, 2, 3]);
        phases = phaseMap;
    end
    
    in.specs = in.specs.*exp(-1i*phaseArray);
    
    if ~exist('weightMap', 'var')
        %initalizing variables for the weighting functions
        weights = zeros([in.sz(in.dims.coils), in.sz(in.dims.x), in.sz(in.dims.y)]);
        
        %getting all the first points at each voxel
        for x = 1:in.sz(in.dims.x)
            for y = 1:in.sz(in.dims.y)
                for i = 1:in.sz(in.dims.coils)
                    %calculate the phase of each each coordinate for all the
                    %coils
                    weights(i,x,y) = in.specs(samplePoint,i,x,y);
                end
            end
        end
        %Get the root sum squared value of coils
        weight_sum = squeeze(rssq(weights, 1));
        %permute for easy multiplication
        weights = permute(weights, [2,3,1]);
        %element wise multiplication
        weights_norm = weights./weight_sum;
        %permute back
        weights_norm = permute(weights_norm, [3,1,2]);
        %getting the weights factor which is first fid
        
        weights = weights_norm;
    else
        weights = weightMap;
    end
    %adding weights to each coil
    permuted_specs = permute(in.specs, [in.dims.coils,in.dims.x,in.dims.y,in.dims.t]);
    combinedSpecs = permuted_specs.*weights;
    combinedSpecs = permute(combinedSpecs, [4,1,2,3]);
    combinedSpecs = squeeze(sum(combinedSpecs, in.dims.coils));
    

    %updating MRSI parameters
    in.dims.x = 2;
    in.dims.y = 3;
    in.dims.coils = 0;
    in.flags.addedrcvrs = 1;
    in.specs = combinedSpecs;
    in.sz = size(in.specs);

end