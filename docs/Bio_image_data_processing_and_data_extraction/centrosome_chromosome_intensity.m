function Centrosome(filename)
% clc;
clearvars -except filename;
close all;

% definition of a guassian function
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';

% read out the data file and make a variable (for both channels)
Tiff = tiffreadVolume(filename);
ChroImage = double(Tiff(:,:,1));
CentImage = double(Tiff(:,:,2));

% determine the number of frames of the file 
T = size(CentImage, 1);

% find the peaks at a particular frame (input the frame number) and find
% the subsequent peaks within the same region of the previous peaks
[~, locs, ws] = findpeaks(double(CentImage(8,:)), 'SortStr','descend');
I = (ws>5.0);
locs = locs(I);
R = locs(1:2); % the first two largest peaks as the initial guess

CentPos = NaN(T,2); % initialization of the position of the centrosome (value=NaN)
S = 20; % the size of the peak search window and ignore the outside of this window

% forward reading of the frames 
for i=8:T
    % read out and normalize the data for timepoint i
    x = double((1:size(CentImage,2))');
    y = double(CentImage(i,:)');
    y = y./max(y); % normalization (the value of y must not exceed 1)
    
    % centrosome 1
    StartPoint = [0.05, R(1), 10, 0.8]; % initial guess for fitting the guassian (a, b, c, d)
    % define the window region boundary (beginning and final boundary)
    x1 = max(1, floor(R(1))-S);
    x2 = min(length(x), floor(R(1))+S);
    
    % the value of x1 should always be less than x2, if not, the value will
    % be NaN
    if (x1<x2)
        f1 = fit(x(x1:x2), y(x1:x2), gaussEqn, 'Start', StartPoint);
    else
        f1.b = NaN;
    end
    
    % centrosome 2
    StartPoint = [0.05, R(2), 10, 0.8];
    x1 = max(1, floor(R(2))-S);
    x2 = min(length(x), floor(R(2))+S);
    if (x1<x2)
        f2 = fit(x(x1:x2), y(x1:x2), gaussEqn, 'Start', StartPoint);
    else
        f2.b = NaN;
    end
    
    % if everything works fine, set R to frame as stated for the initial
    % guess (frame 8)
    if (~isnan(f1.b) && ~isnan(f2.b))
        R = [f1.b, f2.b];
    else
        R = CentPos(8,:);
    end
        
    CentPos(i,:) = [f1.b, f2.b]; % save the centrosome position 
end

% backward reading of the frames 
R = CentPos(8,:);
for i=8:-1:1
    x = double((1:size(CentImage,2))');
    y = double(CentImage(i,:)');
    y = y./max(y);

    StartPoint = [0.05, R(1), 10, 0.8];
    f1 = fit(x(floor(R(1))-S:floor(R(1))+S), y(floor(R(1))-S:floor(R(1))+S), gaussEqn, 'Start', StartPoint);

    StartPoint = [0.05, R(2), 10, 0.8];
    f2 = fit(x(floor(R(2))-S:floor(R(2))+S), y(floor(R(2))-S:floor(R(2))+S), gaussEqn, 'Start', StartPoint);

    R = [f1.b, f2.b];

    CentPos(i,:) = R;
end

% sort the centrosome position; centrosome 1 on the left and centrosome 2
% on the right 
if(CentPos(1,1)>CentPos(1,2))
    J = CentPos(:,1);
    CentPos(:,1) = CentPos(:,2);
    CentPos(:,2) = J;
end


%% Chromosome detection

% the chromosome detetction is the same as the centrosome detection  but
% the detection will be within the boundary of the two detected centrosomes

% M=[11:12];close all;figure;y = double(ChroImage(M,:)');y = y./max(y);plot(y)

ChromPos = NaN(T,3);
MetaTime = -1; % time of anaphase

for i=8:T
    
    % define the centrosome positions and normalize the peak for the
    % chromosome position
    X1 = CentPos(i,1);
    X2 = CentPos(i,2);
    y = double(ChroImage(i,:)');
    y = y./max(y);
    
    % find the peaks for the chromsomes
    [~,locs, ws, ps] = findpeaks(y);
    I = find(ws>4 & locs>X1 & locs<X2);
    ps = ps(I); % magnitutde of the peaks
    locs = locs(I);

    if (i==8)
        [~,j] = max(ps);
        locs = locs(j);
    end

    if (ChromPos(i-1, 3)==1 && length(locs)>=2)
        d1 = abs(ChromPos(i-1, 1) - locs(:));
        if (sum(d1>2)<length(d1))
            j = find(d1>=3);
            locs(j) = [];
        end
    end

    % one chromosome spot
    if (length(locs)==1)
        ChromPos(i, 1:2) = [locs, locs];
        ChromPos(i, 3) = 1;
    end

    % two chromosome spots, first timepoint from metaphase to anaphase
    if (length(locs)>=2 && ChromPos(i-1, 3)==1 && MetaTime<0)
        d1 = (ChromPos(i-1, 1) - locs(:));
        [~, j1] = minpositive(d1);
        if (isempty(j1))
            d1 = abs(ChromPos(i-1, 1) - locs(:));
            [~, j1] = minpositive(d1);
        end
        ChromPos(i, 1) = locs(j1);
        locs(j1) = [];

        d2 = (locs(:) - ChromPos(i-1, 2));
        [~, j2] = minpositive(d2);
        if (isempty(j2))
            d2 = abs(locs(:) - ChromPos(i-1, 2));
            [~, j2] = minpositive(d2);
        end
        ChromPos(i, 2) = locs(j2);

        ChromPos(i, 2) = locs(j2);
        if (ChromPos(i, 2)<ChromPos(i, 1))
            ChromPos(i, 1) = j;
            ChromPos(i, 1) = ChromPos(i, 2);
            ChromPos(i, 2) = j;
        end

        ChromPos(i, 3) = 2;
        MetaTime = i;
    end

    % two chromosome spots, cotinue of anaphase
    if (length(locs)>=2 && isnan(ChromPos(i,3)))

        d1 = (ChromPos(i-1, 1) - locs(:));
        [~, j1] = minpositive(d1);
        if (isempty(j1))
            d1 = abs(ChromPos(i-1, 1) - locs(:));
            [~, j1] = minpositive(d1);
        end
        if (isempty(j1))
            fprintf(sprintf('Something wrong with this file: %s\n', filename));
        else

            ChromPos(i, 1) = locs(j1);
            locs(j1) = [];
        end

        d2 = (locs(:) - ChromPos(i-1, 2));
        [~, j2] = minpositive(d2);
        if (isempty(j2))
            d2 = abs(locs(:) - ChromPos(i-1, 2));
            [~, j2] = minpositive(d2);
        end
        if (isempty(j2))
            fprintf(sprintf('Something wrong with this file: %s\n', filename));
        else
            ChromPos(i, 2) = locs(j2);
        end

        if (~isempty(j1) && ~isempty(j2))
            ChromPos(i, 3) = 2;
            if (d1(j1)>20 || d2(j2)>20)
                ChromPos(i, :) = NaN;
                break;
            end
        end
    end

    if (ChromPos(i, 3)==1 && ChromPos(i-1, 3)==2)
        ChromPos(i, :) = NaN;
        break;
    end
end

for i=8:-1:1
    X1 = CentPos(i,1);
    X2 = CentPos(i,2);
    [~,locs, ws] = findpeaks(double(ChroImage(i,:)));
    I = find(ws>5 & locs>X1 & locs<X2);
    locs = locs(I);

    % one chromosome spot
    if (length(locs)==1)
        ChromPos(i, 1:2) = [locs, locs];
        ChromPos(i, 3) = 1;
    end

    if (length(locs)==2)
        d1 = abs(ChromPos(i+1, 1) - locs(:));
        [~, j1] = minpositive(d1);
        if (isempty(j1))
            fprintf(sprintf('Something wrong with this file: %s\n', filename));
        else
            ChromPos(i, 1:2) = [locs(j1), locs(j1)];
            ChromPos(i, 3) = 1;
        end
    end

    if (abs(ChromPos(i,1)-ChromPos(i+1,1))>20 && isempty(locs))
        ChromPos(i, :) = NaN;
        break;
    end
end

m = find(~isnan(ChromPos(:,1)));
m = m(end);
if(ChromPos(m,1)>ChromPos(m,2))
    J = ChromPos(:,1);
    ChromPos(:,1) = ChromPos(:,2);
    ChromPos(:,2) = J;
end

D = ChromPos(:,2) - ChromPos(:,1);
[~, maxD] = max(D);
lastD = [];
for i=maxD:T
    if (D(i)<0.9*max(D))
        lastD = i;
        break;
    end
end
ChromPos(lastD:T,:) = NaN;

figure('Visible','off');
plot((1:T),CentPos(:,1),'b-', (1:T),CentPos(:,2),'r-', 'LineWidth', 5);
hold on;
plot((1:T),ChromPos(:,1),'k-', (1:T),ChromPos(:,2),'k-', 'LineWidth', 5);
xlabel('time'); ylabel('position');
pbaspect([1 1 1]);
filename = sprintf('%s.png', filename(1:length(filename)-4));
saveas(gcf,filename)

AllPos = [(1:T)', CentPos, ChromPos];
filename = sprintf('%s.csv', filename(1:length(filename)-4));
writematrix(AllPos, filename);

end

function [mins, idxes] = minpositive(Array)
Array(Array<0) = nan;
if (sum(isnan(Array))==length(Array))
    mins = [];
    idxes = [];
else
    [mins, idxes] = min(Array);
end
end

