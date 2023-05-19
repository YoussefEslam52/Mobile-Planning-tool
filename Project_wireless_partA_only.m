%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
%%%annotated project in pdf and final_wireless_project_code%%%
%Constants
channels = 340;
f = 900; % Frequency in MHz *10^9
hb = 20; % BS height in meters
hm = 1.5; % MS height in meters
sensitivity = -95; % MS sensitivity in dBm
traffic_per_user = 0.025; % Traffic intensity per user in Erlangs
path_loss_exponent = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask user for input parameters
GOS = input('GOS (ex:0.02): ');
city_area = input('The City Area (In km^2): ');
user_density = input('User Density (users/km^2): ');
SIRmin_dB = input('Minimum SIR Required (In dB): ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Choosing The Method Of Sectorization
input1 = 1;
while input1
    sectorization_method = input(['Please Enter The Sectorization method: ...' ...
        '\n a = omni-directional \n b = 120° sectorization  \n c = 60° sectorization \n '], 's');
    sectorization_method = lower(sectorization_method); % Convert input to lowercase
    if sectorization_method == 'a'
        sectors = 1;
        i0 = 6;
        input1 = 0;
    elseif sectorization_method == 'b'
        sectors = 3;
        i0 = 2;
        input1 = 0;
    elseif sectorization_method == 'c'
        sectors = 6;
        i0 = 1;
        input1 = 0;
    else
        disp('Wrong parameter')
        input1 = 1;
    end
end
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cluster Size Calculations
SIR_ratio = 10^(SIRmin_dB/10); % SIR ratio in dB
Q = (SIR_ratio*i0)^(1/path_loss_exponent);
A_TEMP = ((Q+1).^2)/3;
Cluster_size = ceil(A_TEMP);
found_n = false;  % flag to indicate if a solution has been found
while ~found_n
    for i = 0:Cluster_size
        for k = 0:Cluster_size
            if i^2 + i*k + k^2 == Cluster_size
                found_n = true;
                break;
            end
        end
        if found_n
            break;
        end
    end
    if ~found_n
        Cluster_size = Cluster_size + 1;
    end
end
disp(['Cluster_size = ', num2str(Cluster_size)]);%->1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Intensity and Number of cells

num_channels_per_sector = floor(channels / (Cluster_size * sectors));

%Solving The Erling B equation using fzero function

fun = @(A) GOS - (A^num_channels_per_sector/factorial(num_channels_per_sector)) ...
/ sum(A.^((0:num_channels_per_sector))./factorial(0:num_channels_per_sector)); 
traffic_intensity_per_sector = fzero(fun, [0, 1000]);

traffic_intensity_per_cell = traffic_intensity_per_sector * sectors;
total_traffic_intensity = user_density * city_area * traffic_per_user;
Num_cells = ceil((total_traffic_intensity) / traffic_intensity_per_cell);
disp(['Total Number of Cells = ', num2str(Num_cells)]);%->2
disp(['Traffic Intensity per Cell = ', num2str(traffic_intensity_per_cell), ' Erlang']);%->4
disp(['Traffic Intensity per Sector = ', num2str(traffic_intensity_per_sector), ' Erlang']);%->5
cell_Area = city_area/Num_cells;
cell_red = sqrt((2*cell_Area)/(3*sqrt(3)));
disp(['Cell radius = ', num2str(cell_red),' Km']);%->3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the transmitted power
L= Hata(f, hm, hb, cell_red);
Ptx = sensitivity+ L ;
disp(['Ptx = ', num2str(Ptx),' dBm']);%->6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the MS received power as a function of distance
D=Q*cell_red;
d = 0:0.01:D;
L= Hata(f, hm, hb, d);
Prx = Ptx - L ;
plot(d, Prx, 'LineWidth', 2);%->7
xlabel('Distance from BS (km)');
ylabel('Received Power (dBm)');
title('Received Power vs. Distance');
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=Hata(f,hm,hb,distance)
CH =0.8+(1.1*log10(f)-0.7)*hm-1.5 *log10(f);
L=69.55 + 26.16 * log10(f) - 13.82 * log10(hb) - CH + (44.9 - 6.55 * log10(hb)) * log10(distance);
end