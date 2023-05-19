%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1)Cluster Size Vs SIR%
% Range of SIRmin_dB values to plot

SIRmin_dB_range = 1:0.001:30;
path_loss_exponent = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize arrays to store cluster sizes for each sectorization method
cluster_sizes_a = zeros(size(SIRmin_dB_range)); %
cluster_sizes_b = zeros(size(SIRmin_dB_range));
cluster_sizes_c = zeros(size(SIRmin_dB_range));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(SIRmin_dB_range)
    SIRmin_dB = SIRmin_dB_range(i);

    % Calculate SIR ratio in dB
    SIR_ratio = 10^(SIRmin_dB/10); 

    % Calculate the cluster size for omni-directional sectorization
    i0 = 6;
    cluster_sizes_a(i) = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent);
    % Calculate the cluster size for 120° sectorization
    i0 = 2 ;
    cluster_sizes_b(i) = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent);
    % Calculate the cluster size for 60° sectorization
    i0 = 1;
    cluster_sizes_c(i) = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the cluster size versus SIRmin_dB for all sectorization methods
figure;
plot(SIRmin_dB_range, cluster_sizes_a, 'LineWidth', 2);
hold on;
plot(SIRmin_dB_range, cluster_sizes_b, 'LineWidth', 2);
plot(SIRmin_dB_range, cluster_sizes_c,  'LineWidth', 2);
ylabel('Cluster Size');
xlabel('Minimum SIR Required (dB)');
title('Cluster Size vs. Minimum SIR Required for Different Sectorization Methods');
legend('Omni-directional', '120° sectorization', '60° sectorization', 'Location', 'Northwest');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(2)&(3)Number of cells and Intensity Vs GOS


SIRmin_dB = input(['Number of cells and Intensity Vs GOS \n...' ...
    ' Minimum SIR Required(In dB): ']); % Use 14dB then use 19dB
%Constants
city_area = 100;%km^2
user_density= 1400;%users/km^2
traffic_per_user = 0.025; % Traffic intensity per user in Erlangs
path_loss_exponent = 4;
channels = 340;
% Range of GOS_values values to plot
GOS_values = 0.01:0.01:0.3;
SIR_ratio = 10^(SIRmin_dB/10); 
% Initialize arrays to store cluster sizes for each sectorization method
n_i0 = [6 2 1]; n_sectors = [1 3 6];
Num_cells = zeros(length(n_i0), length(GOS_values));
traffic_intensity_per_cell = zeros(length(n_i0), length(GOS_values));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for X = 1 : length(n_i0)
    i0 = n_i0(X); sectors = n_sectors(X);
    Cluster_size = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent);
    num_channels_per_sector = floor(channels / (Cluster_size * sectors));
    for Y = 1:length(GOS_values)
        GOS = GOS_values(Y);

        num_channels_per_sector = floor(channels / (Cluster_size * sectors));
    %Solving The Erling B equation using fzero function 
        fun = @(A) GOS - (A^num_channels_per_sector/factorial(num_channels_per_sector)) ...
        / sum(A.^((0:num_channels_per_sector))./factorial(0:num_channels_per_sector)); 
        traffic_intensity_per_sector = fzero(fun, [0, 1000]); 


        traffic_intensity_per_cell(X,Y) = traffic_intensity_per_sector * sectors;
        total_traffic_intensity = user_density * city_area * traffic_per_user;
        Num_cells(X,Y) = ceil((total_traffic_intensity) / traffic_intensity_per_cell(X,Y));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the number of cells vs GOS
figure;
plot(GOS_values, Num_cells(1,:), 'r', GOS_values, Num_cells(2,:),...
    'g', GOS_values, Num_cells(3,:), 'b', 'LineWidth', 2);
title(sprintf('Number of Cells vs GOS (SIRmin = %d dB)', SIRmin_dB));
xlabel('GOS');
ylabel('Number of Cells');
legend('Omni directional', '120° sectorization', '60° sectorization');
grid on;
% Plot the traffic intensity per cell vs GOS
figure;
plot(GOS_values, traffic_intensity_per_cell(1,:), 'r', GOS_values,...
    traffic_intensity_per_cell(2,:), 'g', GOS_values, traffic_intensity_per_cell(3,:), 'b', 'LineWidth', 2);
title(sprintf('Traffic Intensity per Cell vs GOS (SIRmin = %d dB)', SIRmin_dB));
xlabel('GOS');
ylabel('Traffic Intensity per Cell');
legend('Omni directional', '120° sectorization', '60° sectorization');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(4)&(5) Cell Redius and Number of cells Vs User Density


SIRmin_dB = input(['Cell Redius and Number of cells Vs User Density \n ...' ...
    ' Minimum SIR Required(In dB): ']); % Use 14dB then use 19dB
%Constants
city_area = 100;%km^2
GOS=0.02;
traffic_per_user = 0.025; % Traffic intensity per user in Erlangs
path_loss_exponent = 4;
channels = 340;
user_density = 100:2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIR_ratio = 10^(SIRmin_dB/10); 
% Initialize arrays to store cluster sizes for each sectorization method
n_i0 = [6 2 1]; n_sectors = [1 3 6];
Num_cells = zeros(length(n_i0), length(user_density));
cell_red = zeros(length(n_i0), length(user_density));
for X = 1 : length(n_i0)
    i0 = n_i0(X); sectors = n_sectors(X);
    Cluster_size = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent);
    num_channels_per_sector = floor(channels / (Cluster_size * sectors));
%Solving The Erling B equation using fzero function
      fun = @(A) GOS - (A^num_channels_per_sector/factorial(num_channels_per_sector)) ...
        / sum(A.^((0:num_channels_per_sector))./factorial(0:num_channels_per_sector)); 
        traffic_intensity_per_sector = fzero(fun, [0, 1000]); 

        traffic_intensity_per_cell = traffic_intensity_per_sector * sectors;
        for Y = 1:length(user_density)
            total_traffic_intensity = user_density(Y) * city_area * traffic_per_user;
            Num_cells(X,Y) = ceil((total_traffic_intensity)) / traffic_intensity_per_cell;
            cell_Area = city_area/Num_cells(X,Y);
            cell_red(X,Y) = sqrt((2*cell_Area)/(3*sqrt(3)));
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the number of cells versus user density
figure();
plot(user_density, Num_cells(1,:), 'b', 'LineWidth', 2);
hold on;
plot(user_density, Num_cells(2,:), 'r', 'LineWidth', 2);
plot(user_density, Num_cells(3,:), 'g', 'LineWidth', 2);
grid on;
xlabel('User density (users/km^2)', 'FontSize', 12);
ylabel('Number of cells', 'FontSize', 12);
legend('Omni directional', '120° sectorization', '60° sectorization');
title(['Number of cells vs. user density (SIR_{min} = ' num2str(SIRmin_dB) ' dB)'], 'FontSize', 14);
% Plot the cell radius versus user density
figure();
plot(user_density, cell_red(1,:), 'b', 'LineWidth', 2);
hold on;
plot(user_density, cell_red(2,:), 'r', 'LineWidth', 2);
plot(user_density, cell_red(3,:), 'g', 'LineWidth', 2);
grid on;
xlabel('User density (users/km^2)', 'FontSize', 12);
ylabel('Cell radius (km)', 'FontSize', 12);
legend('Omni directional', '120° sectorization', '60° sectorization');
title(['Cell radius vs. user density vs. user density ...' ...
    ' (SIR_{min} = ' num2str(SIRmin_dB) ' dB)'], 'FontSize', 14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=Hata(f,hm,hb,distance)
CH =0.8+(1.1*log10(f)-0.7)*hm-1.5 *log10(f);
L=69.55 + 26.16 * log10(f) - 13.82 * log10(hb) - CH + (44.9 - 6.55 * log10(hb)) * log10(distance);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cluster_size = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent)
    Q = (SIR_ratio*i0)^(1/path_loss_exponent);
    A_TEMP = ((Q+1).^2)/3;
    Cluster_size = ceil(A_TEMP);
    found_n = false;  % flag to indicate if a solution has been found
    while ~found_n
        for j = 0:Cluster_size
            for k = 0:Cluster_size
                if j^2 + j*k + k^2 == Cluster_size
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%