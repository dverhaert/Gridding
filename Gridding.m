%% (C) 2017  Danilo Verhaert/TU Delft
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
load('UVWmatlab');

%% MODE SELECT
% mode = 1 for simple, mode = 2 for w-projection and mode = 3 for
% interpolation.
mode = 1;

% Size of support kernel
oversample = 8; % Number of times we're oversampling
supportsize = 3; % Support for gridding function in pixels
cellsize = 1.08*13107.2/2048; % Cellsize of output grid in wavelengths as defined in Romein. cellsize in Cornwell is 50.
gridsize = 2048; % Size of the grid
nchan = 16; % Number of channels

% Note that these next three variables are set to match with the dataset of
% visibilities. In the basic case, there are a certain number of
% baselines, all measuring in certain timesteps and blocks of these 
% timesteps. When your visibilities are not defined by baselines,
% timesteps and blocks, either set the multiple of these variables to the
% number of points in your dataset or edit the code so it loops from 1
% until your number of data points.

nbaselines = 946; % Number of baselines
ntimesteps = 20; % Number of timesteps present in the visibilities
nblocks = 100; % Number of 'blocks' of visibility data

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

%% Initfrequencies

% Radio frequency spectrum we want to use for UVW projection is around 30 MHz - 10 GHZ
% Romein uses +-60 MHz
% Cornwell uses +-1,4 GHz
for ch = 1:nchan
    % Tweak this frequency depending on your baseline
    frequencies(ch) = (59908828.7353515625 + 12207.03125 *(ch-1))/speed_of_light;
end

%% initVisibilities

% Just for the sake of testing all visibilities are set to 1 at the moment.
% These can of course be loaded from the database.\

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

assert(ntimesteps<=2160);

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
        if(k == wsize/2)
            for j = 1:csize
                j2 = ((j-ccenter)/oversample)^2;
                for i = 1:csize
                    r2 = j2 + ((i-ccenter)/oversample)^2;
                    support(i,j,k) = exp(-r2);
                end
            end
        else
            w = k-wsize/2;
            fScale = sqrt(abs(w)*wcellsize*frequencies(1))/cellsize;
            for j = 1:csize
                j2 = ((j-ccenter)/oversample)^2;
                for i = 1:csize
                    r2 = j2 + ((i-ccenter)/oversample)^2;
                    support(i,j,k) = exp(-r2/(k*fScale));
                end
            end
        end
    end
    
    % Interpolation gridding
elseif mode ==3
    % Generates the csize*csize*wsize kernel support
    % Note that csize is much smaller in comparison to oversampling mode.
    for k = 1:wsize
        if(k == wsize/2)
            for j = 1:csize
                j2 = (j-ccenter)^2;
                for i = 1:csize
                    r2 = j2 + ((i-ccenter))^2;
                    support(i,j,k) = exp(-r2);
                end
            end
        else
            w = k-wsize/2;
            fScale = sqrt(abs(w)*wcellsize*frequencies(1))/cellsize;
            for j = 1:csize
                j2 = ((j-ccenter))^2;
                for i = 1:csize
                    r2 = j2 + ((i-ccenter))^2;
                    support(i,j,k) = exp(-r2/(k*fScale));
                end
            end
        end
    end
    
end


%% Simple mode
if mode == 1
    % These two lines are here so that zeroes can easily be removed
    UVWuvzeros(:,1:2) = UVW(1:nblocks*ntimesteps*nbaselines,1:2);
    UVWuv = UVWuvzeros;
    % Or if you want remove all zero entrys:
    % UVWuv = UVWuvzeros(any(UVWuvzeros,2),:);
    
    % Preallocate memory for speedup
    UVWu_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWv_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    

    
    for ch = 1: nchan
        UVWu_scaled(:,ch) = frequencies(ch)*UVWuv(:,1)/cellsize;
        UVWv_scaled(:,ch) = frequencies(ch)*UVWuv(:,2)/cellsize;
        %axis([-1024, 1024, -1024, 1024]);
        % Plot v against u
        plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
        hold on;
    end
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
                    
                    uscaled = UVWu_scaled(((block-1) * ntimesteps + time-1)*nbaselines + bl, ch);
                    % Matlab command 'fix' rounds to zero.
                    iu = fix(uscaled);
                    fracu = fix(oversample*(uscaled-iu));
                    % Index of the point when rounded and moved to the
                    % center of the grid
                    iu = iu + gridsize/2;
                    % Index required for support matrix
                    uindex = fracu + ccenter;
                    
                    vscaled = UVWv_scaled(((block-1) * ntimesteps + time-1)*nbaselines + bl, ch);
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
    
elseif mode == 2
    UVWuvwzeros(:,:) = UVW(1:nblocks*ntimesteps*nbaselines,:);
    UVWuvw = UVWuvwzeros;
    
    %Or if you want remove all zero entrys:
    %UVWuvw = UVWuvwzeros(any(UVWuvwzeros,2),:);
    
    UVWu_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWv_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;
    UVWw_scaled(nblocks*ntimesteps*nbaselines,nchan) = 0;

    for ch = 1: nchan
        UVWu_scaled(:,ch) = frequencies(ch)*UVWuvw(:,1)/cellsize;
        UVWv_scaled(:,ch) = frequencies(ch)*UVWuvw(:,2)/cellsize;
        UVWw_scaled(:,ch) = frequencies(ch)*UVWuvw(:,3)/wcellsize;
        plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
        hold on;
    end
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
    
    for ch = 1: nchan
        UVWu_scaled(:,ch) = frequencies(ch)*UVWuvw(:,1)/cellsize;
        UVWv_scaled(:,ch) = frequencies(ch)*UVWuvw(:,2)/cellsize;
        UVWw_scaled(:,ch) = frequencies(ch)*UVWuvw(:,3)/wcellsize;
        plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
        hold on;
    end
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
% exact same visibilities, since when there are a lot of zeroes the
% center of your dataset this becomes hard to retrieve except by using
% an advanced algorithm.
% if mode == 1
%     UVWuvzeros(:,1:2) = UVW(1:nblocks*ntimesteps*nbaselines,1:2);
%     UVWuv = UVWuvzeros;
%     
%     %Or if you want remove all zero entries:
%     %UVWuv = UVWuvzeros(any(UVWuvzeros,2),:);
%     
%     for ch = 1: nchan
%         UVWu_scaled(:,ch) = frequencies(ch)*UVWuv(:,1)/cellsize;
%         UVWv_scaled(:,ch) = frequencies(ch)*UVWuv(:,2)/cellsize;
%     end
%     
%     data(1:nblocks*ntimesteps*nbaselines, nchan) = 0;
%     
%     fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
%     fprintf('Degridding...\n');
%     % Over all data points (defined as blocks*baselines*timesteps) and for
%     % each corresponding channel
%     for block = 1:nblocks
%         for bl = 1:nbaselines
%             for time = 1:ntimesteps
%                 for ch = 1:nchan
%                     
%                     sumviswt = 0;
%                     i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
%                     
%                     uscaled = UVWu_scaled(i, ch);
%                     % Matlab command 'fix' rounds to zero.
%                     iu = fix(uscaled);
%                     fracu = fix(oversample*(uscaled-iu));
%                     % Index of the point when rounded and moved to the
%                     % center of the grid
%                     iu = iu + gridsize/2;
%                     % Index required for support matrix
%                     uindex = fracu + ccenter;
%                     
%                     vscaled = UVWv_scaled(i, ch);
%                     iv = fix(vscaled);
%                     fracv = fix(oversample*(vscaled-iv));
%                     iv = iv + gridsize/2;
%                     vindex = fracv + ccenter;
%                     
%                     for suppv = -supportsize:supportsize
%                         for suppu = -supportsize:supportsize
%                             data(i,ch) =  data(i,ch) + support(uindex + suppu,vindex + suppv)*grid(iv + suppv, iu + suppu);
%                             sumviswt = sumviswt + support(uindex + suppu,vindex + suppv);
%                         end
%                     end
%                     data(i,ch) = data(i,ch)/sumviswt;
%                 end
%             end
%         end
%     end
%     
%     %% Degrid a W-projected UV-grid
% elseif mode == 2
%     UVWuvwzeros(:,:) = UVW(1:nblocks*ntimesteps*nbaselines,:);
%     UVWuvw = UVWuvwzeros;
%     
%     %Or if you want remove all zero entrys:
%     %UVWuvw = UVWuvwzeros(any(UVWuvwzeros,2),:);
%     
%     for ch = 1: nchan
%         UVWu_scaled(:,ch) = frequencies(ch)*UVWuvw(:,1)/cellsize;
%         UVWv_scaled(:,ch) = frequencies(ch)*UVWuvw(:,2)/cellsize;
%         UVWw_scaled(:,ch) = frequencies(ch)*UVWuvw(:,3)/wcellsize;
%     end
%     
%     data(1:nblocks*ntimesteps*nbaselines, nchan) = 0;
%     
%     fprintf('The program will go over a total of %d data samples for every channel. \n', nblocks*ntimesteps*nbaselines);
%     fprintf('Degridding...\n');
%     % Over all data points (defined as blocks*baselines*timesteps) and for
%     % each corresponding channel
%     for block = 1:nblocks
%         for bl = 1:nbaselines
%             for time = 1:ntimesteps
%                 for ch = 1:nchan
%                     
%                     sumviswt = 0;
%                     i = ((block-1) * ntimesteps + time-1)*nbaselines + bl;
%                     
%                     uscaled = UVWu_scaled(i, ch);
%                     % Matlab command 'fix' rounds to zero.
%                     iu = fix(uscaled);
%                     fracu = fix(oversample*(uscaled-iu));
%                     % Index of the point when rounded and moved to the
%                     % center of the grid
%                     iu = iu + gridsize/2;
%                     % Index required for support matrix
%                     uindex = fracu + ccenter;
%                     
%                     vscaled = UVWv_scaled(i, ch);
%                     iv = fix(vscaled);
%                     fracv = fix(oversample*(vscaled-iv));
%                     iv = iv + gridsize/2;
%                     vindex = fracv + ccenter;
%                     
%                     wscaled = UVWw_scaled(i, ch);
%                     windex = fix(wscaled) + wsize/2;
%                     
%                     for suppv = -supportsize:supportsize
%                         for suppu = -supportsize:supportsize
%                             data(i,ch) =  data(i,ch) + support(uindex + suppu,vindex + suppv, windex)*grid(iv + suppv, iu + suppu);
%                             sumviswt = sumviswt + support(uindex + suppu,vindex + suppv, windex);
%                         end
%                     end
%                     data(i,ch) = data(i,ch)/sumviswt;
%                 end
%             end
%         end
%     end
% end

%% No degridding for interpolation implemented (yet).
fprintf('All done!\n');

fprintf('Time spent: %d seconds\n', toc);
