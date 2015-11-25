%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Description:    Sea Ice Concentration Anomaly Data analysis:
%                   Task 2: The dataset is split into 3 equal parts of nine
%                   years each and then the correlation based graphs are
%                   constructed. The following subtasks are performed for
%                   each Graph constructed:
%                       1. Plot histogram of Degree distribution
%                       2. Identify Super Nodes
%                       3. Calculate Clustering Coefficient
%                       4. Calculate Characteristic Path Length
%                       5. Compare (3) & (4) for a Random Graph
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc
clear;
addpath(fullfile(pwd, 'matlab_bgl'));   %%% Add library path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dataset Path
dataDir = './Large_Dataset/with_landmask';

%%% R Threshold
R_THRESH = 0.90;
% R_THRESH = 0.95;

%%% Super Nodes XLS File
SuperNodes_XLS_File = ['SuperNodes_Task_2_Yr1_R_', num2str(R_THRESH), '.xls'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Dataset Parameters
%%% [Select Year Range]
YearRange = 1979:1987;
%YearRange = 1988:1996;
%YearRange = 1997:2005;

WeekRange = 1:52;
NumRows = 448;
NumCols = 304;
NumYear = size(YearRange, 2);
NumWeek = size(WeekRange, 2);
DATA_VAL_LAND = 168;
DATA_VAL_MISS = 157;

CORR_BLOCK_LEN = 10000;
WRITE_XLS_ENABLE = true;

%%% Read the Dataset from binary file
i = 1;
Data_Raw_Vect = zeros(NumRows*NumCols, NumYear*NumWeek);
for Year = YearRange
    for Week = WeekRange
        if(Week < 10)
            WeekStr = ['0', num2str(Week)];
        else
            WeekStr = num2str(Week);
        end
        
        FileName = ['/', num2str(Year), '/', 'diff', ...
                    'w', WeekStr, 'y', num2str(Year), '+', 'landmask'];
        FullFilePath = fullfile(dataDir, FileName);
        FileID = fopen(FullFilePath);
        Data_Raw_Vect(:, i) =  fread(FileID, 'float', 'ieee-le');
        i = i + 1;
        fclose(FileID);
    end
end
clearvars Year Week WeekStr FileName FullFilePath FileID i;


%%% Create Land & Missing data MAP and corresponding Mask
Land_Miss_Mask = (Data_Raw_Vect(:,1)==DATA_VAL_LAND) .* DATA_VAL_LAND;
Land_Miss_Mask = Land_Miss_Mask + (Data_Raw_Vect(:,1)==DATA_VAL_MISS) .* DATA_VAL_MISS;

Land_Miss_Map = zeros(NumRows, NumCols, 3);
Map_Temp = logical(reshape(Land_Miss_Mask, [NumCols NumRows])');
Land_Miss_Map(:,:,1) = 255 * Map_Temp;
Land_Miss_Map(:,:,2) = 255 * Map_Temp;
Land_Miss_Map(:,:,3) = 255 * Map_Temp;
clearvars Map_Temp;


%%% Remove unwanted data point from the Raw data
Data_Vect = Data_Raw_Vect(~Land_Miss_Mask, :);
NumDataPoints = size(Data_Vect, 1);
clearvars Data_Raw_Vect;


%%% Construct the Correlation-based Graph using Pearson Correlation
Total_Edges = 0;
%%% Binary Sparse matrix to represent the graph
Corr_Graph = logical( sparse(NumDataPoints, NumDataPoints) );
for x = 1:CORR_BLOCK_LEN:(NumDataPoints-1)
    %%% CORR: Matlab inbuilt function to calculate linear corr coeff
    r = corr(Data_Vect(x:min((x+CORR_BLOCK_LEN-1),end),:)', Data_Vect((x+1):end,:)', 'type', 'Pearson');
    r = abs(r)';
    for i = 1:size(r, 2)
        r_i = find(r(i:end,i) >= R_THRESH) + x + (i-1);
        if( ~isempty(r_i) )
            Corr_Graph(x+(i-1), r_i) = logical(true);
            Total_Edges = Total_Edges + length(r_i);
        end
    end
end
Corr_Graph = Corr_Graph + Corr_Graph';
clearvars r;


%%% Find Degree for all the Vertices
Degree_Temp = sum(Corr_Graph, 2);
if( exist('Degree', 'var') )
    clearvars Degree;
end
[Degree(:,1), ~, Degree(:,2)] = find(Degree_Temp);
NumVert = size(Degree, 1);
clearvars Degree_Temp;


%%% Plot the histogram of Degree Distribution
figure(1);
histogram(Degree(:,2), 'BinMethod', 'integers');
clear title; title('Histogram of Degree Distribution');


%%% Identifying the Super Nodes:
%%% Calc Mean Degree
Degree_Mean = mean(Degree(:,2));

%%% Arrange Nodes in 2D and Find Super nodes
Nodes_Temp = zeros(NumDataPoints, 1);
Nodes_Temp(Degree(:,1), :) = Degree(:,2);
Nodes_2D = zeros(NumRows*NumCols, 1);
Nodes_2D( ~Land_Miss_Mask ) = Nodes_Temp;
Nodes_2D = reshape(Nodes_2D, [NumCols NumRows])';
Super_Nodes_2D = (Nodes_2D>Degree_Mean) .* Nodes_2D;
clearvars Nodes_Temp;

%%% Map the Degree values to colormap 
ColorMap = parula(range(Nodes_2D(:))+1);
ColorMap(1,:) = [0 0 0];
Nodes_RGB = ind2rgb(Nodes_2D - min(Nodes_2D(:)) + 1, ColorMap);
Super_Nodes_RGB = ind2rgb(Super_Nodes_2D - min(Super_Nodes_2D(:)) + 1, ColorMap);
Super_Red_Map = Land_Miss_Map;
Super_Red_Map(:,:,1) = Land_Miss_Map(:,:,1) + (Super_Nodes_2D ~= 0) * 255;

%%% Show the Degree Distribution overlaid on Map of the Earth
figure(2);
image( Super_Red_Map );
clear title; title('Super Nodes: Highlighted in Red');
figure(3);
image( Land_Miss_Map + Nodes_RGB );
clear title; title('All Nodes: Degree Density Representation');
figure(4);
image( Land_Miss_Map + Super_Nodes_RGB );
clear title; title('Super Nodes: Degree Density Representation');

%%% Store the Supernodes in a Spreadsheet
if WRITE_XLS_ENABLE
    [Super_Nodes_XY(:,1), Super_Nodes_XY(:,2), Super_Nodes_XY(:,3)] = find(Super_Nodes_2D);
    try
        if( exist( SuperNodes_XLS_File, 'file') )
            delete( SuperNodes_XLS_File );
        end
        title = {'Super Node (Row)' 'Super Node (Col)' 'Degree'};
        xlswrite(SuperNodes_XLS_File, [title; num2cell(Super_Nodes_XY)]);
    catch
       disp('Cannot write file'); 
    end
end


%%% Calculate Clustering Coeff and Characteristic Path Length 
index = 1;
Total_Pairs = 0;
Short_Path_Sum = 0;
Cluster_E = zeros(size(Degree, 1), 1);
for v=1:NumDataPoints
    N = find(Corr_Graph(v, :))';
    if(isempty(N))
        continue;
    end
    E_v = 0;
    for x=1:size(N, 1)
        for y=x+1:size(N, 1)
            if(Corr_Graph(N(x,1), N(y,1)) == true)
                E_v = E_v + 1;
            end
        end
    end
    %%% e(v)
    Cluster_E(index, 1) = E_v;
    
    %%% Finding Shortest Paths from all nodes to all the other nodes
    %%% **! SHORTEST_PATHS: uses MatlabBGL package by David Gleich !**
    Dist = shortest_paths(double(Corr_Graph), v)';
    Dist( isinf(Dist) ) = [];
    Short_Path_Sum = Short_Path_Sum + sum( Dist );
    Total_Pairs = Total_Pairs + size(Dist, 2) - 1;
    
    index = index + 1;
end
clearvars N E_v;

Cluster_Coeff_v = (2 * Cluster_E(:,1)) ./ (Degree(:,2) .* (Degree(:,2)-1));
Cluster_Coeff_Size = size(Cluster_Coeff_v, 1); 
Cluster_Coeff_v( isnan(Cluster_Coeff_v) ) = [];
Cluster_Coeff_v( isinf(Cluster_Coeff_v) ) = [];

%%% Clustering Coeffcient:
Cluster_Coeff = sum(Cluster_Coeff_v) / Cluster_Coeff_Size;
%%% Characteristic Path Length:
CharPath_Len = Short_Path_Sum / Total_Pairs;

%%% Clustering Coeffcient (Random Graph):
Cluster_Coeff_Rand = Degree_Mean / NumVert;
%%% Characteristic Path Length (Random Graph):
CharPath_Len_Rand = log(NumVert) / log(Degree_Mean);

%%% Comparision:
Cluster_Coeff_Ratio = Cluster_Coeff / Cluster_Coeff_Rand;
CharPath_Len_Ratio = CharPath_Len / CharPath_Len_Rand;


%%% Results:
fprintf('============================================');
fprintf('\nResuts:\n');
fprintf('============================================');
fprintf('\nDegree Mean = %f\n', Degree_Mean);

fprintf('\nClustering Coeff (G_r) =\t %f', Cluster_Coeff);
fprintf('\nClustering Coeff (G_rand) =\t %f', Cluster_Coeff_Rand);
fprintf('\nRatio (G_r / G_rand) =\t\t %f\n', Cluster_Coeff_Ratio);

fprintf('\nCharacteristic Path Length (L_r) =\t %f', CharPath_Len);
fprintf('\nCharacteristic Path Length (L_rand) = %f', CharPath_Len_Rand);
fprintf('\nRatio (L_r / L_rand) =\t\t\t %f', CharPath_Len_Ratio);
fprintf('\n');
