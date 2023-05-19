clear; clc;
plot_cluster_size_vs_SIR()%->1
[Num_cell, traffic_intensity_per_cell] = nocells_Apercell_vs_GOS(14);%->2
[Num_cell, traffic_intensity_per_cell] = nocells_Apercell_vs_GOS(19);%->3
[Num_cell, cell_red] = cellred_nocells_vs_user_dens(14);%->4
[Num_cell, cell_red] = cellred_nocells_vs_user_dens(19);%->5
function [Num_cells, traffic_intensity_per_cell] = nocells_Apercell_vs_GOS(SIRmin_dB)

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
n_i0 = [6 2 1]; 
n_sectors = [1 3 6];
Num_cells = zeros(length(n_i0), length(GOS_values));
traffic_intensity_per_cell = zeros(length(n_i0), length(GOS_values));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for X = 1 : length(n_i0)
    i0 = n_i0(X); sectors = n_sectors(X);
    Cluster_size = calculate_cluster_size(i0, SIR_ratio, path_loss_exponent);
    num_channels_per_sector = floor(channels / (Cluster_size * sectors));
    for Y = 1:length(GOS_values)
        GOS = GOS_values(Y);
    %Solving TheErling B equation using fzero function 
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

end
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
function [Num_cells, cell_red] = cellred_nocells_vs_user_dens(SIRmin_dB)

%Constants
city_area = 100; %km^2
GOS = 0.02;
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

    % Solvingthe Erlang B equation using fzero function
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
title(['Number ofcells vs. user density (SIR_{min} = ' num2str(SIRmin_dB) ' dB)'], 'FontSize', 14);

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
title(['Cell radius vs. user density ...' ...
    ' (SIR_{min} = ' num2str(SIRmin_dB) ' dB)'], 'FontSize', 14);

end
function plot_cluster_size_vs_SIR()

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

end
