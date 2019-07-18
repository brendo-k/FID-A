%io_loadspec_twix.m
%Jamie Near, McGill University 2014.
%Edits from Franck Lamberton, 2017.
%
% USAGE:
% out=io_loadspec_twix(filename);
% 
% DESCRIPTION:
% Reads in siemens twix raw data (.dat file) using the mapVBVD.m and 
% twix_map_obj.m functions from Philipp Ehses (philipp.ehses@tuebingen.mpg.de).
% 
% op_loadspec_twix outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Siemens twix data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadspec_twix(filename);


%read in the data using the new mapVBVD.  This code has been adapted to 
%handle both single RAID files and multi-RAID files.  The vast majority of
%Siemens twix data comes as a single RAID file, but I've encoundered a few 
%multi-RAID files, particularly when using VD13D.  The way to distinguish
%them here is that a for a single RAID file, mapVBVD will output a struct, 
%whereas for a multi-RAID file, mapVBVD will output a cell array of structs.
%This code assumes that the data of interest is in the last element of the 
%cell array (possibly a bad assumption under some circumstances):
twix_obj=mapVBVD(filename);
if isstruct(twix_obj)
    disp('single RAID file detected.');
    RaidLength=1;
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    twix_obj=twix_obj{RaidLength};
end
dOut.data=twix_obj.image();
version=twix_obj.image.softwareVersion;
sqzSize=twix_obj.image.sqzSize; 
sqzDims=twix_obj.image.sqzDims;


%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging 
%to find the correct sequence name.  
sequence=twix_obj.hdr.Config.SequenceFileName;  

%Try to find out what sequnece this is:
isSpecial=~isempty(strfind(sequence,'rm_special')) ||...  %Is this Ralf Mekle's SPECIAL sequence?
            ~isempty(strfind(sequence,'vq_special'));  %or the CIBM SPECIAL sequence?
isjnSpecial=~isempty(strfind(sequence,'jn_svs_special')) ||...  %or Jamie Near's SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Adiab_Special')); %or Masoumeh Dehghani's Adiabatic SPECIAL sequence?
isjnMP=~isempty(strfind(sequence,'jn_MEGA_GABA')); %Is this Jamie Near's MEGA-PRESS sequence?
isjnseq=~isempty(strfind(sequence,'jn_')) ||... %Is this another one of Jamie Near's sequences 
        ~isempty(strfind(sequence,'md_'));      %or a sequence derived from Jamie Near's sequences (by Masoumeh Dehghani)?
isWIP529=~isempty(strfind(sequence,'edit_529')); %Is this WIP 529 (MEGA-PRESS)?
isWIP859=~isempty(strfind(sequence,'edit_859')); %Is this WIP 859 (MEGA-PRESS)?
isMinn=~isempty(strfind(sequence,'eja_svs_')); %Is this one of Eddie Auerbach's (CMRR, U Minnesota) sequences?
isSiemens=(~isempty(strfind(sequence,'svs_se')) ||... %Is this the Siemens svs PRESS seqeunce?
            ~isempty(strfind(sequence,'csi_se')) ||... %Or the Siemens CSI PRESS sequence?
            ~isempty(strfind(sequence,'svs_st')) || ... % Or the Siemens svs STEAM sequence?
            ~isempty(strfind(sequence,'csi_st'))) && ... % Or the Siemens CSI STEAM sequence?
            isempty(strfind(sequence,'eja_svs'));    %And make sure it's not 'eja_svs_steam', which also contains the string 'svs_st'.
        
%Is this a CSI sequence?:
isCSI=~isempty(strfind(sequence,'csi')) ||...
        ~isempty(strfind(sequence,'CSI'));

%If this is the SPECIAL sequence, it probably contains both inversion-on
%and inversion-off subspectra on a single dimension, unless it is the VB
%version of Jamie Near's SPECIAL sequence, in which case the subspecs are
%already stored on separate dimensions.  
%Both Ralf Mekle's SPECIAL and the VD-VE version of Jamie Near's SPECIAL sequence 
%do not store the subspectra along a separate dimension of the data array, 
%so we will separate them artifically:
%25 Oct 2018: Due to a recent change, the VE version of Jamie Near's MEGA-PRESS 
%sequence also falls into this category. 
if isSpecial ||... %Catches Ralf Mekle's and CIBM version of the SPECIAL sequence 
        (strcmp(version,'vd') && isjnSpecial) ||... %and the VD/VE versions of Jamie Near's SPECIAL sequence
        (strcmp(version,'vd') && isjnMP);  %and the VD/VE versions of Jamie Near's MEGA-PRESS sequence                                                   
    squeezedData=squeeze(dOut.data);
    if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
        data(:,:,:,1)=squeezedData(:,:,[1:2:end-1]);
        data(:,:,:,2)=squeezedData(:,:,[2:2:end]);
        sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/2 2];
    elseif twix_obj.NCol>1 && twixObj.image.NCha==1
        data(:,:,1)=squeezedData(:,[1:2:end-1]);
        data(:,:,2)=squeezedData(:,[2:2:end]);
        sqzSize=[sqzSize(1) sqzSize(2)/2 2];
    end
    if isjnseq
        sqzDims{end+1}='Set';
    else
        sqzDims{end+1}='Ida';
    end
else
    data=dOut.data;
end

%Squeeze the data to remove singleton dims
fids=squeeze(data);

%noticed that in the Siemens PRESS and STEAM sequences, there is sometimes
%an extra dimension containing unwanted reference scans or something.  Remove them here.
if isSiemens && (strcmp(version,'vd') || strcmp(version,'ve')) && strcmp(sqzDims{end},'Phs')
    sqzDims=sqzDims(1:end-1);
    sqzSize=sqzSize(1:end-1);
    if ndims(fids)==4
        fids=fids(:,:,:,2);
        fids=squeeze(fids);
    elseif ndims(fids)==3
        fids=fids(:,:,2);
        fids=squeeze(fids);
    elseif ndims(fids)==2
        fids=fids(:,2);
        fids=squeeze(fids);
    end
end

%Make a pulse sequence identifier for the header (out.seq);
seq=sequence;

%Find the magnetic field strength:
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;

%Find the number of averages:
Naverages=twix_obj.hdr.Meas.Averages;

%Find out if multiple coil elements were used:
Ncoils=twix_obj.hdr.Meas.iMaxNoOfRxChannels;  

%Find the TE:
TE = twix_obj.hdr.MeasYaps.alTE{1};  %Franck Lamberton

%Find the TR:
TR = twix_obj.hdr.MeasYaps.alTR{1};  %Franck Lamberton

%Now begin indexing the dimensions of the data array. ie. create the dims
%structure, which specifies which dimensions of the data array are being
%used to hold the time-domain data, the multiple coil channels, the
%average, the sub-spectra, and any additional dimensions.
sqzDims_update=sqzDims;
dimsToIndex=[1:length(sqzDims)];

%First index the dimension of the time-domain data
dims.t=find(strcmp(sqzDims,'Col'));
if ~isempty(dims.t)
    %remove the time dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.t);
else
    dims.t=0;
    error('ERROR:  Spectrom contains no time domain information!!');
end

%Now index the dimension of the coil channels
dims.coils=find(strcmp(sqzDims,'Cha'));
if ~isempty(dims.coils)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.coils);
else
    dims.coils=0;
end

%Now index the dimension of the averages
if strcmp(version,'vd') || strcmp(version,'ve')
    if isMinn
        dims.averages=find(strcmp(sqzDims,'Set'));
    else
        dims.averages=find(strcmp(sqzDims,'Ave'));
    end
else
    dims.averages=find(strcmp(sqzDims,'Set'));
end
if ~isempty(dims.averages)
    %remove the averages dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
else
    %If no Averages dimension was found, then check for a "Repetitions"
    %dimension.  If that is found, store it under "averages".  If both
    %"Averages" and "Repetitions" dimensions are found, "Repetitions" will
    %be indexed under "Extras", since "Repetitions is not currently an
    %option in FID-A.
    dims.averages=find(strcmp(sqzDims,'Rep'));
    if ~isempty(dims.averages)
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
    else
        %If neither an "Averages" or a "Repetitions" dimension is found,
        %then set the FID-A "Averages" dimension to zero.
        dims.averages=0;
    end
end

%Now index the CSI dimensions
dims.csi_x=find(strcmp(sqzDims,'Lin'));
dims.csi_y=find(strcmp(sqzDims,'Phs'));
if ~isempty(dims.csi_x)
    %remove the csi_x dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.csi_x);
else 
    %If a "csi_x" dimension is not found,
    %then set the FID-A "csi_x" dimension to zero.
    dims.csi_x=0;
end
if ~isempty(dims.csi_y)
    %remove the csi_y dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.csi_y);
else 
    %If a "csi_y" dimension is not found,
    %then set the FID-A "csi_y" dimension to zero.
    dims.csi_y=0;
end


%Now we have indexed the dimensions containing the timepoints, the coil
%channels, the averages, and the x- and y- CSI encoding.  As we indexed 
%each dimension, we removed the corresponding index from the dimsToIndex 
%vector.  At this point, if there are any values left in the dimsToIndex 
%vector, then there must be some additional dimensions that need indexing.  
%We assume that if sub-spectra exist, then these must be indexed in either 
%the 'Ida' dimension (for all Jamie Near's VB-version pulse sequences), the 
%'Set' dimension (for all Jamie Near's VD/VE-version pulse sequences), the 
%'Eco' dimension (for the WIP 529 MEGA-PRESS sequence or the Minnesota 
%MEGA-PRESS sequence), or the 'Ide' dimension (for the WIP 859 MEGA-PRESS 
%sequence). 
if ~isempty(dimsToIndex)
    %Now index the dimension of the sub-spectra
    if isjnseq  || isSpecial
        if strcmp(version,'vd') || strcmp(version,'ve')
            dims.subSpecs=find(strcmp(sqzDims,'Set'));
        else
            dims.subSpecs=find(strcmp(sqzDims,'Ida'));
        end
    elseif isWIP529 || isMinn
        dims.subSpecs=find(strcmp(sqzDims,'Eco'));
    elseif isWIP859
        dims.subSpecs=find(strcmp(sqzDims,'Ide'));
    else
        dims.subSpecs=dimsToIndex(1);
    end
    if ~isempty(dims.subSpecs)
        %remove the sub-spectra dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.subSpecs);
    else
        dims.subSpecs=0;
    end
else
    dims.subSpecs=0;
end

%And if any further dimensions exist after indexing the sub-spectra, call
%these the 'extras' dimension.  
if ~isempty(dimsToIndex)
    %Now index the 'extras' dimension
    dims.extras=dimsToIndex(1);
    if ~isempty(dims.extras)
        %remove the extras dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.extras);
    else
        dims.extras=0;
    end
else
    dims.extras=0;
end

%Now that we've indexed the dimensions of the data array, we now need to
%permute it so that the order of the dimensions is standardized:  we want
%the order to be as follows:  
%   1) time domain data.  
%   2) coils.
%   3) averages.
%   4) subSpecs.
%   5) csi_x.
%   6) csi_y.
%   7) extras.
if length(sqzDims)==7
    fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.csi_x dims.csi_y dims.extras]);
    dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=5;dims.csi_y=6;dims.extras=7;
elseif length(sqzDims)==6
    if dims.extras==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=5;dims.csi_y=6;dims.extras=0;
    elseif dims.csi_y==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.csi_x dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=5;dims.csi_y=0;dims.extras=6;
    elseif dims.csi_x==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.csi_y dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=0;dims.csi_y=5;dims.extras=6;
    elseif dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=4;dims.csi_y=5;dims.extras=6;
    elseif dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=2;dims;averages=0;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=5;dims.extras=6;
    elseif dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=5;dims.extras=6;
    end
elseif length(sqzDims)==5
    if dims.extras==0 && dims.csi_y==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subspecs dims.csi_x]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=5;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subspecs dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=0;dims.csi_y=5;dims.extras=0;
    elseif dims.extras==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=4;dims.csi_y=5;dims.extras=0;
    elseif dims.extras==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=5;dims.extras=0;
    elseif dims.extras==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=5;dims.extras=0;
    elseif dims.csi_y==0 && dims.csi_x==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subspecs dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=0;dims.csi_y=0;dims.extras=5;
    elseif dims.csi_y==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.csi_x dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=4;dims.csi_y=0;dims.extras=5;
    elseif dims.csi_y==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.csi_x dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=0;dims.extras=5;
    elseif dims.csi_y==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.csi_x dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=0;dims.extras=5;
    elseif dims.csi_x==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.csi_y dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=4;dims.extras=5;
    elseif dims.csi_x==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.csi_y dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=0;dims.csi_y=4;dims.extras=5;
    elseif dims.csi_x==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.csi_y dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=0;dims.csi_y=4;dims.extras=5;
    elseif dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=3;dims.csi_y=4;dims.extras=5;
    elseif dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.csi_x=3;dims.csi_y=4;dims.extras=5;
    elseif dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.csi_x=3;dims.csi_y=4;dims.extras=5;
    end
elseif length(sqzDims)==4
    if dims.extras==0 && dims.csi_y==0 && dims.csi_x==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.csi_x]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=4;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.csi_x]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.csi_x]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=4;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=2;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=4;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=0;dims.csi_y=4;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=0;dims.csi_y=4;dims.extras=0;
    elseif dims.extras==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=3;dims.csi_y=4;dims.extras=0;
    elseif dims.extras==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.csi_x=3;dims.csi_y=4;dims.extras=0;
    elseif dims.extras==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.csi_x=3; dims.csi_y=4;dims.extras=0;
    elseif dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=0; dims.csi_y=0;dims.extras=4;
    elseif dims.csi_y==0 && dims.csi_x==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=0; dims.csi_y=0;dims.extras=4;
    elseif dims.csi_y==0 && dims.csi_x==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=0; dims.csi_y=0;dims.extras=4;
    elseif dims.csi_y==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.csi_x dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=3; dims.csi_y=0;dims.extras=4;
    elseif dims.csi_y==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.csi_x dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.csi_x=3; dims.csi_y=0;dims.extras=4;
    elseif dims.csi_y==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs dims.csi_x dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.csi_x=3; dims.csi_y=0;dims.extras=4;
    elseif dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.csi_x dims.csi_y dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=2; dims.csi_y=3;dims.extras=4;
    end
elseif length(sqzDims)==3
    if dims.extras==0 && dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.csi_x==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.csi_x==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.csi_x]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=3;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.csi_x]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=3;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs dims.csi_x]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.csi_x=3;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.csi_y]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=3;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=3;dims.extras=0;
    elseif dims.extras==0 && dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.csi_x dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=2;dims.csi_y=3;dims.extras=0;
    elseif dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=0;dims.extras=3;
    elseif dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=0;dims.extras=3;
    elseif dims.csi_y==0 && dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.csi_x dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=2;dims.csi_y=0;dims.extras=3;
    elseif dims.csi_x==0 && dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.csi_y dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=2;dims.extras=3;
    end
elseif length(sqzDims)==2
    if dims.extras==0 && dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.csi_x==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.csi_x=0;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_y==0 && dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.csi_x]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=2;dims.csi_y=0;dims.extras=0;
    elseif dims.extras==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.csi_y]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=2;dims.extras=0;
    elseif dims.csi_y==0 && dims.csi_x==0 && dims.subSpecs==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.csi_x=0;dims.csi_y=0;dims.extras=2;
    end
end

%Now get the size of the data array:
sz=size(fids);

%Now take fft of time domain to get specs:
specs=fftshift(ifft(fids,[],dims.t),dims.t);

%Now take fft of csi dimensions to get image data.  NOTE:  Before this step
%We should make sure that the centre of k-space corresponds to the centre
%of the csi_x and csi_y dimensions.  If not, we will have to shift the data
% prior to this FT step.  Fix at later date.  Also note: K-space will not be
%stored in the FID-A structure.  If needed, we'll have to transform back to
%k-space.;
if dims.csi_x
    fids=ifft(fids,[],dims.csi_x);
    specs=ifft(specs,[],dims.csi_x);
    kx=[];  %Fill this out later
end
if dims.csi_y
    fids=ifft(fids,[],dims.csi_y);
    specs=ifft(specs,[],dims.csi_y);
    ky=[];  %Fill this out later
end
    

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
spectralwidth=1/dwelltime;
    
%Get TxFrq
txfrq=twix_obj.hdr.Meas.Frequency;

%Get Date
%date = getfield(regexp(twix_obj.hdr.MeasYaps.tReferenceImage0, ...
%'^".*\.(?<DATE>\d{8})\d*"$', 'names'), 'DATE');  %Franck Lamberton

date=''; %The above code for extracting the date from the header 
         %was causing problems.  Since date is not critical
         %for almost any applications, removing it now to be fixed at a
         %later date.

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
    else
        averages=sz(dims.subSpecs);
        rawAverages=1;
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
    else
        averages=1;
        rawAverages=1;
    end
end

%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

%Find the number of points acquired before the echo so that this
%information can be stored in the .pointsToLeftshfit field of the data
%structure.  Depending on the pulse sequence used to acquire the data, the
%header location of this parameter is different.  For product PRESS
%seqeunces, the value is located in twix_obj.image.freeParam(1).  For WIP
%sequences, the value is located in twix_obj.image.cutOff(1,1).  For CMRR
%sequences, the value is located in twix_obj.image.iceParam(5,1).  Special
%thanks to Georg Oeltzschner for decoding all of this and sharing the
%information with me:

if isWIP529 || isWIP859
    leftshift = twix_obj.image.cutOff(1,1);
elseif isSiemens
    leftshift = twix_obj.image.freeParam(1);
elseif isMinn
    leftshift = twix_obj.image.iceParam(5,1);
else
    leftshift = twix_obj.image.freeParam(1);
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.65;

t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=seq;
out.te=TE/1000;
out.tr=TR/1000;
out.pointsToLeftshift=leftshift;


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end
out.flags.isCSI=isCSI;



%DONE
