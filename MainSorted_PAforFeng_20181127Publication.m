%% Main
clear
clc

FolderThisAnalysis = 'S:\HCS_Platform\Data\PaulAntony\Feng\MitoGraph\BrunoSantos\20181127Publication';
mkdir(FolderThisAnalysis)
f_LogDependencies(mfilename, FolderThisAnalysis)

fileDetails = dir('S:\HCS_Platform\Data\BrunoSantos\20180703 new stainings n1\*.czi')
fileDetails = sortrows(struct2table(fileDetails), 'bytes')
PathCell = cell(height(fileDetails), 1);
for i = 1:height(fileDetails)
    fileDetailThis = fileDetails(i,:);
    PathCell{i} = [fileDetailThis.folder{:}, filesep, fileDetailThis.name{:}];
end

fileDetails = [fileDetails, PathCell];
fileDetails.Properties.VariableNames{end} = 'Path';
files = fileDetails.Path;

ObjectsAll = {};

for i = 1:size(files, 1)
    
    try
        %% load images
        fileThis = files{i};
        ImCell = bfopen(fileThis);
        ImCellClipped = ImCell{1,1};
        ImTable = cell2table(ImCellClipped);
        planeBioformats = regexp(ImTable{:, 'ImCellClipped2'}, '.*plane (.{1,3})\/.*;', 'tokens');
        planeBioformats = cellfun(@(x) str2double(x{:}{:}),  planeBioformats, 'UniformOutput', false);
        ImTable = [ImTable, planeBioformats];
        ImTable = sortrows(ImTable, 'Var3');

        ch1_MTR = [];
        plane = 0;
        for p = [1:4:height(ImTable)]
            plane = plane + 1;
            ch1_MTR(:,:,plane) = ImTable{p, 'ImCellClipped1'}{:};
        end
        % vol(ch1_MTR)

        ch2_Tom20 = [];
        plane = 0;
        for p = [2:4:height(ImTable)]
            plane = plane + 1;
            ch2_Tom20(:,:,plane) = ImTable{p, 'ImCellClipped1'}{:};
        end
        % vol(ch2_Tom20)

        ch3_H = [];
        plane = 0;
        for p = [3:4:height(ImTable)]
            plane = plane + 1;
            ch3_H(:,:,plane) = ImTable{p, 'ImCellClipped1'}{:};
        end
        % vol(ch3_H)

        ch4_Tuj1 = [];
        plane = 0;
        for p = [4:4:height(ImTable)]
            plane = plane + 1;
            ch4_Tuj1(:,:,plane) = ImTable{p, 'ImCellClipped1'}{:};
        end
        % vol(ch4_Tuj1)

        %% image analysis
        mkdir([FolderThisAnalysis, filesep, 'Previews'])
        save([FolderThisAnalysis, '\files.mat'], 'files');
        Objects = ImageAnalysisBrunoFeng_20181127Publication(ch1_MTR, ch2_Tom20, ch3_H, ch4_Tuj1, fileThis, FolderThisAnalysis, i);
        ObjectsAll{i} = Objects;
        
        %% Checkpointing
        save([FolderThisAnalysis, filesep, 'data_', num2str(i), '.mat'], 'Objects');
        
        catch e
            Errors{i} = e;
        end

end

Data = vertcat(ObjectsAll{:});
DataCompact = Data(:, ~ismember(Data.Properties.VariableNames, 'NodeDegreeVector'));
save([FolderThisAnalysis, filesep, 'data.mat']);
writetable(DataCompact,  [FolderThisAnalysis, filesep, 'data.xlsx']);
writetable(DataCompact,  [FolderThisAnalysis, filesep, 'data.csv']);

