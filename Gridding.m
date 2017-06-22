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
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Initialize and reset

close all
clear all
% load('UVW1blockmatlab');
load('UVWmatlab');

%% Instead of using 1 kernel, we use 8x8xw_planes kernels, explaining the
% three extra dimensions in the support variable. In this case we proceeed
% first by extending the algorithm to contain w-projection.
%% MODE SELECT
% mode = 1 for simple, mode = 2 for oversample, mode = 3 for interpolate
mode = 1;
% gridding = 0 for degridding, gridding = 1 for gridding
gridding = 1;

% Variables used:
% Gridtype, support, visibilities, UVW, frequencies, supportPixelsUsed. All
% declared later.

GRID_U = 2048;
GRID_V = 2048;

% Size of kernel
SUPPORT_U = 16;
SUPPORT_V = 16;
X = SUPPORT_U;
SUPHALF_U = SUPPORT_U/2;
SUPHALF_V = SUPPORT_V/2;

W_PLANES = 2;

% Oversample subgrid by 8
OVERSAMPLE_U = 8;
OVERSAMPLE_V = 8;

CELL_SIZE_U = 1.08*13107.2/GRID_U;
CELL_SIZE_V = 1.08*13107.2/GRID_V;
CELL_SIZE_W = 8192/W_PLANES;

NR_STATIONS = 44;
BASELINES = (NR_STATIONS * (NR_STATIONS - 1) / 2); % 946

CHANNELS = 16;
TIMESTEPS = 20;
BLOCKS = 100; % Our short data version has 1 blocks (20k rows), was 108, 100

SPEED_OF_LIGHT = 299792458;

% Testvariable kip
kip = 2;

%% Initsupport
% Simple mode
if mode == 1
    % Generates the 16x16 kernel support
    for v = 1:SUPPORT_V
        for u = 1:SUPPORT_U
            support(v,u) = min(v,SUPPORT_V-(v-1))*min(u,SUPPORT_U-(u-1));
        end
    end
    
    % Oversampling mode
    % Generates the 16*16*32*8*8 kernel support
elseif mode == 2
    for w = 1:W_PLANES
        for ou = 1:OVERSAMPLE_U
            for ov = 1:OVERSAMPLE_V
                for v = 1:SUPPORT_V
                    for u = 1:SUPPORT_U
                        support(w,ov,ou,v,u) = w*min(v,SUPPORT_V-(v-1))*min(u,SUPPORT_U-(u-1)) + (w*ov*ou)*1i ;
                    end
                end
            end
        end
    end
    
    % Interpolate mode
    % Generates the 16*16*32 kernel support
elseif mode == 3
    for w = 1:W_PLANES
        for v = 2:SUPPORT_V+1
            for u = 2:SUPPORT_U+1
                % Same as simple but then *w and /normalization
                % *3/(8*17*17)
                support(w,v,u) = (w-1)*min(v-1,SUPPORT_V-((v-1)-1))*min(u-1,SUPPORT_U-((u-1)-1))/(W_PLANES*(SUPPORT_V+1)*(SUPPORT_U+1)/4);
            end
            support(w,v,1) = 0;
            support(w,v,SUPPORT_U+2) = 0;
        end
        for u = 1:SUPPORT_U+2
            support(w,1,u) = 0;
            support(w,SUPPORT_V+2,u) = 0;
        end
    end
end


%% initfrequencies

for ch = 1:CHANNELS
    frequencies(ch) = 59908828.7353515625 + 12207.03125 * (ch-1);
end


% TODO: For degridding, we do not want to run initVisibilities. But for
% now, it's ok.

%% initVisibilities

vis = 2 + i;

for bl = 1:BASELINES
    for time = 1:TIMESTEPS
        for ch = 1:CHANNELS
            visibilities(bl,time,ch) = vis;
        end
    end
end

assert(NR_STATIONS<=44);
assert(TIMESTEPS<=2160);

%% initUVW

for ch = 1:CHANNELS
    scale_u(ch) = frequencies(ch)/(CELL_SIZE_U*SPEED_OF_LIGHT);
    scale_v(ch) = frequencies(ch)/(CELL_SIZE_V*SPEED_OF_LIGHT);
    scale_w(ch) = frequencies(ch)/(CELL_SIZE_W*SPEED_OF_LIGHT);
end

%% TODO: doCPUdegridding

if gridding == 0
    fprintf('Degridding selected. Not implemented yet.');
    
    
    %% doCPUgridding
    
elseif gridding == 1
    % Simple mode
    
    if mode == 1
        UVWuvzeros(:,1:2) = UVW(1:BLOCKS*TIMESTEPS*BASELINES,1:2);
        UVWuv = UVWuvzeros;
        
        %If you want to plot the entire dataset or remove zeroes:
        %UVWuvzeros = UVW(:,1:2);
        %UVWuv = UVWuvzeros(any(UVWuvzeros,2),:);
        
        
        for ch = 1: CHANNELS
            UVWu_scaled(:,ch) = round(scale_u(ch)*UVWuv(:,1));
            UVWv_scaled(:,ch) = round(scale_v(ch)*UVWuv(:,2));
            %axis([-1024, 1024, -1024, 1024]);
            % Plot v against u
            plot(UVWu_scaled(:,ch),UVWv_scaled(:,ch),'.')
            hold on;
        end
        figure;
        
        grid(1:2048,1:2048) = 0;
        count =0;
        
        fprintf('Assigning grid points...\n');
        for block = 1:BLOCKS
            for bl = 1:BASELINES
                for time = 1:TIMESTEPS
                    for ch = 1: CHANNELS
                        % Old grid accessing
                        % grid_u = round(uvw_u(bl,time,ch));
                        % grid_v = round(uvw_v(bl,time,ch));
                        
                        % weight = support(v,u);
                        % % TODO: Figure out localgrid
                        % prod = visibilities(bl,time,ch) * weight;
                        
                        grid_u = UVWu_scaled(((block-1) * TIMESTEPS + time-1)*946 + bl, ch) + GRID_U/2;
                        grid_v = UVWv_scaled(((block-1) * TIMESTEPS + time-1)*946 + bl, ch)+ GRID_V/2;
                        for v = -SUPHALF_V:SUPHALF_V-1
                            for u = -SUPHALF_U:SUPHALF_U-1
                                grid(grid_v+v, grid_u + u) = grid(grid_v+v, grid_u + u) + support(v+SUPHALF_V+1,u+SUPHALF_U+1);
                                
                                % Original CPP code:
                                % for (unsigned pol = 0; pol < POLARIZATIONS; pol ++) {
                                % float2 prod   = visibilities[bl][0][ch][pol] * weight;
                                % float2 &pixel = (*localGrid)[grid_v + v][grid_u + u][pol];
                                % pixel += prod;
                            end
                        end
                    end
                end
            end
        end
        gridflip = flip(grid);
        fprintf('Flipping and plotting...\n');
        %     Y = ifft2(grid);
        %     imagesc(abs(fftshift(Y)));
        %     figure;
        imagesc(gridflip);
        % v is inverted in the grid, flip for clarity
        %surf(grid);
        
        % Oversampling mode
    elseif mode == 2
        for bl = 1:BASELINES
            % TODO: u_end and v_end again set to supportPixelsUsed, or X. Extremely
            % redundant.
            u_end = X;
            v_end = X;
            for v = 1:v_end
                for time = 1:TIMESTEPS
                    for ch = 1: CHANNELS
                        uc = uvw_u(bl,time,ch);
                        vc = uvw_v(bl,time,ch);
                        wc = uvw_w(bl,time,ch);
                        unsignedwc = typecast(uint16(wc),'int16');
                        grid_u = round(uvw_u(bl,time,ch));
                        grid_v = round(uvw_v(bl,time,ch));
                        u_frac = uc - grid_u;
                        v_frac = vc - grid_v;
                        % TODO: Round or floor?
                        ou = OVERSAMPLE_U * u_frac;
                        ov = OVERSAMPLE_V * v_frac;
                        unsignedou = typecast(uint16(ou),'int16');
                        unsignedov = typecast(uint16(ov),'int16');
                        
                        % TODO: Add polarizations. Needs to be added in a lot
                        % of locations.
                        vis = visibilities(bl,time,ch);
                        
                        for u = 1:u_end
                            % Note: The following step takes quite some time
                            % since in this case support is a 5D double!
                            weight = support(unsignedwc+1,unsignedov+1,unsignedou+1,v,u);
                            % TODO: Figure out localgrid
                            % prod = visibilities(bl,time,ch) * weight;
                            % pixel =
                            % pixel = pixel + prod;
                        end
                    end
                end
            end
        end
        
        % Interpolation mode
    elseif mode == 3
        for block = 1:BLOCKS
            for bl = 1:BASELINES
                for time = 1:TIMESTEPS
                    for ch = 1:CHANNELS
                        % This code combines the uvw and currentUVW variables of the
                        % original code
                        uvw_u(bl,time,ch) = scale_u(ch) * UVW(((block-1) * TIMESTEPS + time-1)*946 + bl,1) + GRID_U/2 - 16/2;
                        uvw_v(bl,time,ch) = scale_v(ch) * UVW(((block-1) * TIMESTEPS + time-1)*946 + bl,2) + GRID_V/2 - 16/2;
                        uvw_w(bl,time,ch) = scale_w(ch) * UVW(((block-1) * TIMESTEPS + time-1)*946 + bl,3) + W_PLANES/2;
                        
                        %Testing
                        
                        
                    end
                end
                
            end
        end
        
        for bl = 1:BASELINES-18 % Random number - to test for different grid_u_fracs icm u0 and u1
            for v = 1:X
                scale_u = SUPPORT_U / X;
                scale_v = SUPPORT_V / X;
                for time = 1:TIMESTEPS
                    for ch = 1:CHANNELS
                        % For every data sample,
                        grid_u = uvw_u(bl,time,ch);
                        grid_v = uvw_v(bl,time,ch);
                        w = uvw_w(bl,time,ch);
                        grid_u_int = round(grid_u); % 113,3 --> 113
                        grid_u_frac = grid_u - grid_u_int; % Deze is dan 0,3
                        grid_v_int = round(grid_v);
                        grid_v_frac = grid_v - grid_v_int;
                        
                        support_v = 0.5 + scale_v * (v-grid_v_frac + 0.5);
                        support_v_int = floor(support_v);
                        v1 = support_v - support_v_int;
                        v0 = 1 - v1;
                        support_w_int = floor(w);
                        w1 = w - support_w_int;
                        w0 = 1-w1;
                        
                        vw00 = v0*w0;
                        for u = 1:X
                            support_u = 0.5 + scale_u * (u-grid_u_frac + 0.5);
                            support_u_int = floor(support_u);
                            u1 = support_u - support_u_int;
                            u0 = 1 - u1;
                                % TODO: Finish weighting
%                             weight = weight = u0 * (vw00 * support[support_w_int    ][support_v_int    ][support_u_int    ] +
%                             vw01 * support[support_w_int    ][support_v_int + 1][support_u_int    ] +
%                             vw10 * support[support_w_int + 1][support_v_int    ][support_u_int    ] +
%                             vw11 * support[support_w_int + 1][support_v_int + 1][support_u_int    ]) +
%                             u1 * (vw00 * support[support_w_int    ][support_v_int    ][support_u_int + 1] +
%                             vw01 * support[support_w_int    ][support_v_int + 1][support_u_int + 1] +
%                             vw10 * support[support_w_int + 1][support_v_int    ][support_u_int + 1] +
%                             vw11 * support[support_w_int + 1][support_v_int + 1][support_u_int + 1]);
                            % grid(grid_v+v, grid_u + u) = grid(grid_v+v, grid_u + u) + support(v+SUPHALF_V+1,u+SUPHALF_U+1);
                            
                            
                            % Bij positieve input 0.3 grid_u_frac:
                            % u0 = grid_u_frac, u1 = 1-grid_u_frac
                            % In dit geval dus 0.3 en 0.7.
                            % Bij input -0.3 voor grid_u_frac:
                            % u0 = 1 + grid_u_frac en u1 = -grid_u_frac
                            % In dit geval dus 0.7 en 0.3.
                            % Bij 0 is u0 = 1. Bij 0.01 u0 = 0.01.
                            
                            
                            
                        end
                    end
                end
            end
        end
    end
end

%% printVisibilities (degridding only)

% if gridding == 0
%     printvisibilities = 1;
%
%
%     %% printGrid
%
% elseif gridding == 1
%     count_v = 0;
%     for v = 1:GRID_V
%         count_u = 0;
%         for u = 0:GRID_U
%             % Temporarily commented out until we fount out how grid is
%             % filled
%             %if (grid[v][u][0].x != 0 || grid[v][u][0].y != 0)
%             % We willen grid(v,u) gevuld hebben. Dus 2048*2048.
%             %             if (gridu(v,u) ~= 0 || gridv(v,u) ~= 0)
%             %                 count_u = count_u+1;
%             %                 if(count_u == 0)
%             %                     count_v = count_v+1;
%             %                 end
%             %                 if (count_u<5 && count_v<5)
%             %                     fprintf('CPU: ( %d , %d ): %d \n', u,v,grid(v,u));
%             %                 end
%             %             end
%         end
%
%
%
%     end
% end
%
fprintf('All done!\n');
