function Master(FOLDER)

clc;
close all;
clearvars -except FOLDER;

files = dir(sprintf('%s/*.tif', FOLDER));
size(files)

parfor i=1:size(files,1)
    filename = sprintf('%s/%s', files(i).folder, files(i).name);
    filename
    Centrosome(filename);
end

end