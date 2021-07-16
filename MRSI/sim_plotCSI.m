% sim_plotCSI.m
% Jamie Near, McGill University 2019.
% 
% USAGE:
% sim_plotCSI(in,coilNum)
% sim_plotCSI(in)
% 
% DESCRIPTION:
% This function takes a processed MRSI twix file and plots the data using 
% the matlab plot function. Horizontal coordinates are the x coordinates
% and vertical are y coordinates. The plots at specific (x,y) positions are
% the spectral coordintaes
%
% INPUTS:
% in          = input cell array of simulated spectra from a spatially resolved simulation
% coilNum     = coil to plot if coils have not been combined (default is 1)
%               Use coil 1 if coils have been already combined
%                
%
% OUTPUTS: 
% displays graph of processed MRSI data in twix struct. 


function sim_plotCSI(in, ppmmin, ppmmax, complex_plot, coilNum)

%Argument checks

%check to see if coilNum has been defined
if ~exist('coilNum','var')
    coilNum = 1;
elseif in.dims.coils == 0 && coilNum > 1
    error('no coils to plot')
end

if ~exist('complex_plot','var')
    complex_plot = 'real';
end

if ~exist('ppmmin', 'var')
    if ~isfield(in, 'ppm')
        ppmmin = min(in.t);
    else
        ppmmin = min(in.ppm);
    end    
end

if ~exist('ppmmax', 'var')
    if ~isfield(in, 'ppm')
        ppmmax = max(in.t);
    else
        ppmmax = max(in.ppm);
    end
end

if ~strcmpi(complex_plot, 'real') && ~strcmpi(complex_plot, 'imag') && ~strcmpi(complex_plot, 'abs')
    error("please input 'real', 'imag' or 'abs'")
end

%check if ppm exists to plot
if ~isfield(in, 'ppm')
    range_bool = in.t > ppmmin & in.t < ppmmax;
    ppm = in.t(range_bool);
else
    %if ppm doesn't exist plot time domain
    range_bool = in.ppm > ppmmin & in.ppm < ppmmax;
    ppm = in.ppm(range_bool);
end

specs = permute(in.specs, nonzeros([in.dims.t, in.dims.x, in.dims.y, in.dims.coils, ...
                           in.dims.averages]));
specs = specs(range_bool, :,:);

%get lowercase
complex_plot = lower(complex_plot);

if(strcmp(complex_plot, 'real'))
    yrange = max(max(real(specs),[],1) - min(real(specs),[],1), [], 'all');
elseif (strcmp(complex_plot, 'imag'))
    yrange = max(max(imag(specs),[],1) - min(imag(specs),[],1), [], 'all');
else
    yrange = max(max(abs(specs),[],1) - min(abs(specs),[],1), [], 'all');
end

%scale factors to fit at each (x,y) coordinates
scalefactorX=(0.8*in.deltaX)/(ppm(end)-ppm(1));
scalefactorY=(0.8*in.deltaY)/yrange;
specs = specs.*scalefactorY;
%scale the intensity of the specs to the scalefactorY

%flip the specs along the y axis 

tempSpec = specs;

%create figure and hold 
figure;
hold on;

%set the x and y limits
xlim([(in.xCoordinates(1)-in.deltaX), (in.xCoordinates(end)+in.deltaX)]);
xticks(in.xCoordinates);
ylim([in.yCoordinates(1)-in.deltaX, in.yCoordinates(end)+in.deltaY]);
yticks(in.yCoordinates);

for x = 1:size(specs,in.dims.x)
    for y = 1:size(specs,in.dims.y)
        %first scale the ppm scale so that the range is correct;
        ppm_plot=(ppm-ppm(1))*scalefactorX;
        %now shift it to the correct x-position;
        ppm_plot = ppm_plot + (x-1)*in.deltaX - (0.8*in.deltaX)/2 + in.xCoordinates(1);
        %Now start plotting
        if(strcmp(complex_plot, 'real'))
            plot(ppm_plot, flip(real(tempSpec(:,x,y,coilNum)) - real(tempSpec(1,x,y,coilNum)) + (y-1)*in.deltaY + in.yCoordinates(1)));
        elseif(strcmp(complex_plot, 'imag'))
            plot(ppm_plot, flip(imag(tempSpec(:,x,y,coilNum)) + (y-1)*in.deltaY + in.yCoordinates(1)));
        else
            plot(ppm_plot, flip(abs(tempSpec(:,x,y,coilNum)) + (y-1)*in.deltaY + in.yCoordinates(1)));
        end
        
    end
end
hold off;
end
