function NeurocytoManifold(io_path)
    % NeurocytoManifold: Processes neurocytoskeleton tree data and applies 
    % smoothing algorithms to generate a more accurate 3D model.

    %% User settings and path setup
    parameter_file = [io_path, 'mesh_parameter.txt'];
    input_file = [io_path, 'skeleton_initial.swc'];
    smooth_file = [io_path, 'skeleton_smooth.swc'];

    % Read skeleton information
    neuroTree{1} = load_tree(input_file);
    
    % Load parameters for smoothing
    params = LoadParameter(parameter_file);
    n_noisesmooth = params(1);      % set iteration steps for noise smooth
    ratio_bifur_node = params(2);   % set bifurcation nodes smooth ratio
    ratio_noisesmooth = params(3);  % set noise smooth ratio
    seg_length = params(4);         % set bezier smooth segments length 
    
    % Extract coordinates and connection information
    [sect_point, bif_pt, n_sect, n_bif] = ExtractTreeInfo(neuroTree{1});
    
    %% Noise Smoothing
    neuroTree{1} = NoiseSmoothNeurocyto(neuroTree{1}, bif_pt, sect_point, n_sect, n_bif, n_noisesmooth, ratio_bifur_node, ratio_noisesmooth);

    %% Bezier Smooth
    neuroTree{2} = BezierSmoothNeurocyto(neuroTree{1}, sect_point, n_sect, bif_pt, seg_length);

    %% Output and Visualization
    swc_tree(neuroTree{2}, smooth_file);
    figure(1); clf; xplore_tree(neuroTree{1});
    figure(2); clf; xplore_tree(neuroTree{2});
end

function neuroTree = NoiseSmoothNeurocyto(neuroTree, bif_pt, sect_point, n_sect, n_bif, n_noisesmooth, ratio_bifur_node, ratio_noisesmooth)
    % Noise smoothing logic
    for index_smooth = 1:n_noisesmooth
        % Optimize bifurcation point positions
        for ii = 1:n_bif
            pt_b = [neuroTree.X(bif_pt{ii}(1)), neuroTree.Y(bif_pt{ii}(1)), neuroTree.Z(bif_pt{ii}(1))];  
            pt_i = [neuroTree.X(bif_pt{ii}(2)), neuroTree.Y(bif_pt{ii}(2)), neuroTree.Z(bif_pt{ii}(2))];  
            pt_j = [neuroTree.X(bif_pt{ii}(3)), neuroTree.Y(bif_pt{ii}(3)), neuroTree.Z(bif_pt{ii}(3))];  
            pt_k = [neuroTree.X(bif_pt{ii}(4)), neuroTree.Y(bif_pt{ii}(4)), neuroTree.Z(bif_pt{ii}(4))];  

            % Adjust bifurcation point positions
            pt_b_after = pt_b * (1 - ratio_bifur_node) + (pt_i + pt_j + pt_k) / 3 * ratio_bifur_node;
            neuroTree.X(bif_pt{ii}(1)) = pt_b_after(1);
            neuroTree.Y(bif_pt{ii}(1)) = pt_b_after(2);
            neuroTree.Z(bif_pt{ii}(1)) = pt_b_after(3);
        end
        
        % Smooth noise points along the sections
        for i = 1:n_sect
            points = [neuroTree.X(sect_point{i}), neuroTree.Y(sect_point{i}), neuroTree.Z(sect_point{i}), neuroTree.D(sect_point{i})];
            smoothed_points = NoiseSmooth(points, ratio_noisesmooth, 1);
            neuroTree.X(sect_point{i}) = smoothed_points(:, 1);
            neuroTree.Y(sect_point{i}) = smoothed_points(:, 2);
            neuroTree.Z(sect_point{i}) = smoothed_points(:, 3);
            neuroTree.D(sect_point{i}) = smoothed_points(:, 4);
        end
    end
end

function neuroTree = BezierSmoothNeurocyto(neuroTree, sect_point, n_sect, bif_pt, seg_length)
    % Bezier smoothing logic
    for i = 1:n_sect
        % Process each segment and apply Bezier smoothing
        start_index = sect_point{i}(1);
        end_index = sect_point{i}(end);
        mode = DetermineSmoothingMode(neuroTree, start_index, end_index);
        
        % Bezier smoothing of the branch
        [smoothed_points, tmp_tangent] = BsplineSmooth([neuroTree.X(sect_point{i}), neuroTree.Y(sect_point{i}), neuroTree.Z(sect_point{i})], neuroTree.D(sect_point{i}), seg_length, mode);
        
        % Insert new smoothed points into the tree structure
        neuroTree = InsertSmoothPoints(neuroTree, smoothed_points, tmp_tangent, start_index, end_index);
    end
end

function mode = DetermineSmoothingMode(neuroTree, start_index, end_index)
    % Determine the mode for Bezier smoothing based on the start and end points
    if B_tree(neuroTree)(end_index) && B_tree(neuroTree)(start_index)
        mode = 1;
    elseif T_tree(neuroTree)(end_index)
        mode = 2;
    elseif start_index == 1
        mode = 3;
    else
        mode = 4;
    end
end

function OutputResults(neuroTree, smooth_file)
    % Save the smoothed tree structure and tangents
    swc_tree(neuroTree, smooth_file);
    disp(['Smoothed tree saved to ', smooth_file]);
end
