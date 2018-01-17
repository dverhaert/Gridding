%% (C) 2018  Danilo Verhaert/TU Delft
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

%% Initialize and reset

close all
clear all
% Some test datasets are included in the repository.
load('uvw'); % Dataset name will be uvw
% load('poslocal'); % Dataset name will be poslocal, so add UVW = poslocal
% load('ska'); % Dataset name will be uvw
% load('testvla'); % Dataset name will be uvw
UVW = uvw;

%% MODE SELECT
% mode = 1 for simple, mode = 2 for w-projection and mode = 3 for
% interpolation.
mode = 2;

% Size of support kernel
oversample = 8; % Number of times we're oversampling
supportsize = 3; % Support for gridding function in pixels
gridsize = 2048; % Size of the grid
cellsize = 1.08*13107.2/2048; % Cellsize of output grid in wavelengths as defined in Romein. cellsize in Cornwell is 50.
nchan = 16; % Number of channels

% Note that these next three variables are set to match with the dataset of
% visibilities. In the basic case, there are a certain number of
% baselines, all measuring in certain timesteps and blocks of these
% timesteps. When your visibilities are not defined by baselines,
% timesteps and blocks, you can use the 3 lines below:
nbaselines = 1;
ntimesteps = 1;
nblocks = length(UVW);

% Otherwise, set these:
% Version where we want to split baselines/timesteps/blocks
% nbaselines = 946; % Number of baselines
% ntimesteps = 20; % Number of timesteps present in the visibilities
% nblocks = 170; % Number of 'blocks' of visibility data

speed_of_light = 299792458;
tic; % Start timing

%% Mode specific variables
if mode == 1
    fprintf('Simple mode selected.\n');
elseif mode == 2
    fprintf('W-projection mode selected.\n');
    wsize = 64; % Number of W-planes (only used in W-projection)
    wcellsize = 8192/wsize; % Corresponding cell size
elseif mode == 3
    fprintf('Interpolation mode selected.\n');
    oversample = 1; % Since we're interpolating, we're not oversampling.
    wsize = 64; % Number of W-planes (only used in W-projection)
    wcellsize = 8192/wsize; % Corresponding cell size
end

csize = 2*(supportsize + 1)*oversample+1; % Size of support matrix
support(csize,csize) = 0;
ccenter = (csize+1)/2; % Center of support matrix

%% Initlabda

% Radio frequency spectrum we want to use for UVW projection is around 30 MHz - 10 GHZ
% Romein uses +-60 MHz
% Cornwell uses +-1,4 GHz
for ch = 1:nchan
    % Tweak this frequency depending on your baseline
    labda(ch) = (59908828.7353515625 + 12207.03125 *(ch-1))/speed_of_light;
end

%% initScaling
% Scaling correction to make grid fit. Scaling is necessary since
% grid points are rounded towards the nearest integer, so there have to be
% some significant differences between visibilities to differentiate them
% from eachother.
% For now we scale so that we match to a grid of [-1024,1024]. Basically,
% in doing this we disregard the cellsize/frequency. Comment these two
% lines if you want to set a custom grid size/compare multiple datasets:
scalecorrection = labda(16)*max(max(abs(UVW(:,1)),abs(UVW(:,2))))/cellsize;
correctionfactor = 1000/scalecorrection;

% And in that case, also uncomment this one:
% correctionfactor = 1;

if mode ~= 1
    % Case where all w-coordinates in dataset are 0, except for one or two
    % outliers. In this case, we don't need to use w-projection, so use mode 1. 
    if mean(abs(UVW(:,3))) < 1e-10
        error('No w-coordinates in data. Please use mode 1 for narrow-field of view');
    else
        % wcorrectionfactor = 1;
        wscalecorrection = labda(16)*max(abs(UVW(:,3)))/wcellsize;
        wcorrectionfactor = 31/wscalecorrection; % Correction factor for w
        % Since we want w to be in the range of 1-64
        % This way all number will be inbetween -31 and 31
        % and adding 32 leads to all numbers being inbetween 1 and 63
    end
end

%% initVisibilities

% Just for the sake of testing all visibilities are set to 1 at the moment.
% These can of course be loaded from the database.

visibilities(nbaselines,ntimesteps,nchan) = 0;
vis = 1;

for bl = 1:nbaselines
    for time = 1:ntimesteps
        for ch = 1:nchan
            visibilities(bl,time,ch) = vis;
            %= visibilities((time-1)*nbaselines + bl, ch);
        end
    end
end

%% Initsupport
% Simple gridding
if mode == 1
    % Generates the csize*csize kernel support
    for i = 1:csize
        i2 = ((i-ccenter)/oversample)^2;
        for j = 1:csize
            r2 = i2 + ((j-ccenter)/oversample)^2;
            support(i,j) = exp(-r2);
        end
    end
    
% W-projection gridding
elseif mode == 2
    % Generates the csize*csize*wsize kernel support
    for k = 1:wsize
        w = k-wsize/2;
        % This if statement makes sure that there will be no w-plane that
        % is zero (by removing abs(w) from the equation)
        if(k == wsize/2)
            fScale = sqrt(wcellsize*labda(1))/cellsize;
        else
            fScale = sqrt(abs(w)*wcellsize*labda(1))/cellsize;
        end
        for j = 1:csize
            j2 = ((j-ccenter)/oversample)^2;
            for i = 1:csize
                r2 = j2 + ((i-ccenter)/oversample)^2;
                support(i,j,k) = exp(-r2/(k*fScale));
            end
        end
    end
    
% Interpolation gridding
elseif mode == 3
    % Generates the csize*csize*wsize kernel support
    % Note that csize is much smaller in comparison to oversampling mode.
    for k = 1:wsize
        w = k-wsize/2;
        % This if statement makes sure that there will be no w-plane that
        % is zero (by removing abs(w) from the equation)
        if(k == wsize/2)
            fScale = sqrt(wcellsize*labda(1))/cellsize;
        else
            fScale = sqrt(abs(w)*wcellsize*labda(1))/cellsize;
        end
        for j = 1:csize
            j2 = ((j-ccenter))^2;
            for i = 1:csize
                r2 = j2 + ((i-ccenter))^2;
                support(i,j,k) = exp(-r2/(k*fScale));
            end
        end
    end
end


%% Simple mode
if mode == 1
    % This line is necessary for when wanting to not use all data, or 
    % to split the data up for possible extra functionality
    UVWuv(:,1:2) = UVW(1:nblocks*ntimesteps*nbaselines,1:2);
    
    % Preallocate memory for speedup
    UVWu_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWv_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    % Scale UVW according to labda, and possibly with the correctionfactor
    % if you want it to fit nicely 
    for ch = 1:nchan
        UVWu_scaled(:,ch) = labda(ch)*UVWuv(:,1)/cellsize*correctionfactor;
        UVWv_scaled(:,ch) = labda(ch)*UVWuv(:,2)/cellsize*correctionfactor;
        % Plot v against u
        plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
        hold on;
    end
    axis([-1024, 1024, -1024, 1024]); % Define range of plot
    figure;
        
    grid(1:gridsize,1:gridsize) = 0; % Create a grid 
    
    fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
    fprintf('Assigning grid points...\n');
    % Over all data points (defined as blocks*baselines*timesteps) and for
    % each corresponding channel
    for block = 1:nblocks
        for bl = 1:nbaselines
            for time = 1:ntimesteps
                for ch = 1:nchan
                    
                    % Get the current u coordinate by defining the index iUVWv_scaled(i, ch);
                    i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
                    uscaled = UVWu_scaled(i, ch);
                    
                    % Matlab command 'fix' rounds to zero.
                    iu = fix(uscaled);
                    fracu = fix(oversample*(uscaled-iu));
                    % Index of the point when rounded and moved to the
                    % center of the grid
                    iu = iu + gridsize/2;
                    % Index required for support matrix
                    uindex = fracu + ccenter;
                    
                    vscaled = UVWv_scaled(i, ch);
                    iv = fix(vscaled);
                    fracv = fix(oversample*(vscaled-iv));
                    iv = iv + gridsize/2;
                    vindex = fracv + ccenter;
                    
                    for suppv = -supportsize:supportsize
                        for suppu = -supportsize:supportsize
                            grid(iv + suppv, iu + suppu) =  grid(iv + suppv, iu + suppu) + support(uindex + suppu,vindex + suppv)*visibilities(bl,time,ch);
                        end
                    end
                end
            end
        end
    end
    gridflip = flip(grid);
    fprintf('Flipping and plotting...\n');
    
    imagesc(log(gridflip));
    
%% Oversampling mode 
% Note: Lines which are used multiple times throughout the code are not 
% explained again. Please check mode 1 for this.
elseif mode == 2
    UVWuvw(:,:) = UVW(1:nblocks*ntimesteps*nbaselines,:);
    
    % Preallocate memory for speedup    
    UVWu_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWv_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWw_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    
    for ch = 1: nchan
        UVWu_scaled(:,ch) = labda(ch)*UVWuvw(:,1)/cellsize*correctionfactor;
        UVWv_scaled(:,ch) = labda(ch)*UVWuvw(:,2)/cellsize*correctionfactor;
        UVWw_scaled(:,ch) = labda(ch)*UVWuvw(:,3)/wcellsize*wcorrectionfactor;
        plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
        hold on;
    end
    axis([-1024, 1024, -1024, 1024]);
    figure;
    
    grid(1:gridsize,1:gridsize) = 0;
    
    fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
    fprintf('Assigning grid points...\n');
    % Over all data points (defined as blocks*baselines*timesteps) and for
    % each corresponding channel
    for block = 1:nblocks
        for bl = 1:nbaselines
            for time = 1:ntimesteps
                for ch = 1:nchan                    
                    i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
                    uscaled = UVWu_scaled(i, ch);
                    % Matlab command 'fix' rounds to zero.
                    iu = fix(uscaled);
                    fracu = fix(oversample*(uscaled-iu));
                    % Index of the point when rounded and moved to the
                    % center of the grid
                    iu = iu + gridsize/2;
                    % Index required for support matrix
                    uindex = fracu + ccenter;
                    
                    vscaled = UVWv_scaled(i, ch);
                    iv = fix(vscaled);
                    fracv = fix(oversample*(vscaled-iv));
                    iv = iv + gridsize/2;
                    vindex = fracv + ccenter;
                    
                    wscaled = UVWw_scaled(i, ch);
                    windex = fix(wscaled) + wsize/2;
                    
                    for suppv = -supportsize:supportsize
                        for suppu = -supportsize:supportsize
                            grid(iv + suppv, iu + suppu) =  grid(iv + suppv, iu + suppu) + support(uindex + suppu,vindex + suppv, windex)*visibilities(bl,time,ch);
                        end
                    end
                end
            end
        end
    end
    gridflip = flip(grid);
    fprintf('Flipping and plotting...\n');
    
    imagesc(log(gridflip));
    
%% Interpolation
elseif mode == 3
    UVWuvwzeros(:,:) = UVW(1:nblocks*ntimesteps*nbaselines,:);
    UVWuvw = UVWuvwzeros;
    
    %Or if you want remove all zero entrys:
    %UVWuvw = UVWuvwzeros(any(UVWuvwzeros,2),:);
    UVWu_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWv_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWw_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
        
    
    factor = 64/max(UVWuvw(:,3));
    
    for ch = 1: nchan
        UVWu_scaled(:,ch) = labda(ch)*UVWuvw(:,1)/cellsize*correctionfactor;
        UVWv_scaled(:,ch) = labda(ch)*UVWuvw(:,2)/cellsize*correctionfactor;
        UVWw_scaled(:,ch) = labda(ch)*UVWuvw(:,3)/wcellsize*factor;

        
        plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
        hold on;
    end
    axis([-1024, 1024, -1024, 1024]);
    figure;
    
    grid(1:gridsize,1:gridsize) = 0;
    
    fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
    fprintf('Assigning grid points...\n');
    % Over all data points (defined as blocks*baselines*timesteps) and for
    % each corresponding channel
    for block = 1:nblocks
        for bl = 1:nbaselines
            for time = 1:ntimesteps
                for ch = 1:nchan
                    
                    i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
                    uscaled = UVWu_scaled(i, ch);
                    
                    iu = fix(uscaled);
                    % Determine how much a visibility should contribute to
                    % this grid point or the next one.
                    % u0 gives the factor the current point should be
                    % weighted with, and u1 the next point.
                    ufloor = floor(uscaled);
                    u1 = uscaled - ufloor;
                    u0 = 1 - u1;
                    % Index of the point when rounded and moved to the
                    % center of the grid
                    iu = iu + gridsize/2;
                    % Index required for support matrix
                    uindex = ccenter;
                    
                    vscaled = UVWv_scaled(i, ch);
                    iv = fix(vscaled);
                    % Do the same for v
                    vfloor = floor(vscaled);
                    v1 = vscaled - vfloor;
                    v0 = 1 - v1;
                    iv = iv + gridsize/2;
                    vindex = ccenter;
                    
                    wscaled = UVWw_scaled(i, ch);
                    iw = fix(wscaled);
                    windex = iw + wsize/2;
                    % Do the same for w
                    wfloor = floor(wscaled);
                    w1 = wscaled - wfloor;
                    w0 = 1 - w1;
                    
                    for suppv = -supportsize:supportsize
                        for suppu = -supportsize:supportsize
                            weight = u0*v0*w0*support(uindex + suppu,vindex + suppv, windex) + ...
                                u0*v0*w1*support(uindex + suppu,vindex + suppv, windex + 1) + ...
                                u0*v1*w0*support(uindex + suppu,vindex + suppv + 1, windex) + ...
                                u0*v1*w1*support(uindex + suppu,vindex + suppv + 1, windex + 1) + ...
                                u1*v0*w0*support(uindex + suppu + 1,vindex + suppv, windex) + ...
                                u1*v0*w1*support(uindex + suppu + 1,vindex + suppv, windex + 1) + ...
                                u1*v1*w0*support(uindex + suppu + 1,vindex + suppv + 1, windex) + ...
                                u1*v1*w1*support(uindex + suppu + 1,vindex + suppv + 1, windex + 1);
                            grid(iv + suppv, iu + suppu) =  grid(iv + suppv, iu + suppu) + weight*visibilities(bl,time,ch);
                        end
                    end
                end
            end
        end
    end
    gridflip = flip(grid);
    fprintf('Flipping and plotting...\n');
    
    imagesc(log(gridflip));
    
end

%% Degridding

% Degrid a simple UV-grid. Note that this doesn't give you back the
% exact same visibilities. This because the only input you have is the
% grid, and it is impossible to reconstruct all the visibilities from this,
% only an approximation except by using an advanced algorithm.
degridding = 1;
if degridding == 1;
    if mode == 1
        UVWuv(:,1:2) = UVW(1:length(UVW),1:2);
        
        %Or if you want remove all zero entries:
        %UVWuv = UVWuvzeros(any(UVWuvzeros,2),:);
        
        for ch = 1: nchan
            UVWu_scaled(:,ch) = labda(ch)*UVWuv(:,1)/cellsize*correctionfactor;
            UVWv_scaled(:,ch) = labda(ch)*UVWuv(:,2)/cellsize*correctionfactor;
        end
        
        data(1:nblocks*ntimesteps*nbaselines, nchan) = 0;
        
        fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
        fprintf('Degridding...\n');
        % Over all data points (defined as blocks*baselines*timesteps) and for
        % each corresponding channel
        for block = 1:nblocks
            for bl = 1:nbaselines
                for time = 1:ntimesteps
                    for ch = 1:nchan
                        
                        sumviswt = 0;
                        i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
                        
                        uscaled = UVWu_scaled(i, ch);
                        % Matlab command 'fix' rounds to zero.
                        iu = fix(uscaled);
                        fracu = fix(oversample*(uscaled-iu));
                        % Index of the point when rounded and moved to the
                        % center of the grid
                        iu = iu + gridsize/2;
                        % Index required for support matrix
                        uindex = fracu + ccenter;
                        
                        vscaled = UVWv_scaled(i, ch);
                        iv = fix(vscaled);
                        fracv = fix(oversample*(vscaled-iv));
                        iv = iv + gridsize/2;
                        vindex = fracv + ccenter;
                        
                        for suppv = -supportsize:supportsize
                            for suppu = -supportsize:supportsize
                                data(i,ch) =  data(i,ch) + support(uindex + suppu,vindex + suppv)*grid(iv + suppv, iu + suppu);
                                sumviswt = sumviswt + support(uindex + suppu,vindex + suppv);
                            end
                        end
                        data(i,ch) = data(i,ch)/sumviswt;
                    end
                end
            end
        end
        
%% Degrid a W-projected UV-grid
    elseif mode == 2
        UVWuvw = UVW(1:nblocks*ntimesteps*nbaselines,:);
                
        for ch = 1: nchan
            UVWu_scaled(:,ch) = labda(ch)*UVWuvw(:,1)/cellsize*correctionfactor;
            UVWv_scaled(:,ch) = labda(ch)*UVWuvw(:,2)/cellsize*correctionfactor;
            UVWw_scaled(:,ch) = labda(ch)*UVWuvw(:,3)/wcellsize*wcorrectionfactor;
        end
        
        data(1:nblocks*ntimesteps*nbaselines, nchan) = 0;
        
        fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
        fprintf('Degridding...\n');
        % Over all data points (defined as blocks*baselines*timesteps) and for
        % each corresponding channel
        for block = 1:nblocks
            for bl = 1:nbaselines
                for time = 1:ntimesteps
                    for ch = 1:nchan
                        
                        sumviswt = 0;
                        i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
                        
                        uscaled = UVWu_scaled(i, ch);
                        % Matlab command 'fix' rounds to zero.
                        iu = fix(uscaled);
                        fracu = fix(oversample*(uscaled-iu));
                        % Index of the point when rounded and moved to the
                        % center of the grid
                        iu = iu + gridsize/2;
                        % Index required for support matrix
                        uindex = fracu + ccenter;
                        
                        vscaled = UVWv_scaled(i, ch);
                        iv = fix(vscaled);
                        fracv = fix(oversample*(vscaled-iv));
                        iv = iv + gridsize/2;
                        vindex = fracv + ccenter;
                        
                        wscaled = UVWw_scaled(i, ch);
                        windex = fix(wscaled) + wsize/2;
                        
                        for suppv = -supportsize:supportsize
                            for suppu = -supportsize:supportsize
                                data(i,ch) =  data(i,ch) + support(uindex + suppu,vindex + suppv, windex)*grid(iv + suppv, iu + suppu);
                                sumviswt = sumviswt + support(uindex + suppu,vindex + suppv, windex);
                            end
                        end
                        data(i,ch) = data(i,ch)/sumviswt;
                    end
                end
            end
        end
    end
end
%% No degridding for interpolation (mode 3) implemented (yet).
fprintf('All done!\n');

fprintf('Time spent: %d seconds\n', toc);
