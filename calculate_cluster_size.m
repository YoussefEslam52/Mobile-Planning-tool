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