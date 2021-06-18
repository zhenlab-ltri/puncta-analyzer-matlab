% punctaanalyser Copyright (C) 20071220 Taizo Kawano <tkawano at mshri.on.ca>
%
% This program is free software; you can redistribute it and/modify it 
% under the term of the GNU General Public License as published bythe Free Software Foundation;
% either version 2, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this file.  If not, write to the Free Software Foundation,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

% This program is written for fluoresent puncta analysis of C. elegans dorsal nerve cord.
% Image data prepared by preprocessor, other program written for preprocessing, or Preprocessor_2 ImageJ plugin, are readable data format.
% The simplest usase is; punctaanalyse(`wt-1`) and so on.
% also,  [pcount, plineardensity, pdist, pwidths, pgaps, pintens, pvolumes, fixwidth, fixvol, fixgap] = punctaanalyse(`wt-1`)
% will return all return_values.
% You can set second argument, thresholdcoeff, which involves pcount, plineardensity, pdist, pwidths, pgaps, pintens, pvolumes,
% and third argument, fixthreshold, which involves fixwidth, fixvol, fixgap.
% The default values are 20 and 1000 each. Increasing there result in less puncta ditection.
% The 4th and 5th arguments determin mode (default 0: noisy, 1: silent(batch mode)) and objective (100 or 63).
% Defaultobjective is changed to 0. It will use pixel as out put unit. 091026

function [pcount, plineardensity, pdist, pwidths, pgaps, pintens, pvolumes, fixwidth, fixvol, fixgap, fixwidthlineardensity, fixgaplineardensity]   = punctaanalyser (samplename, thresholdcoeff, fixthreshold, mode, objective)
  % default values.
  defaultthreshold = 20;
  defaultfixthreshold = 1000;
  defaultmode = 0;
  % If you want to change default objective, just change below value, such as defaultobjective = 63;
  %defaultobjective = 100;
  defaultobjective = 0;
  % scalefactor means length of 1 pixel (um). 9.637313 um / 100 pixel with x100 objective, and 15.30746.
  % divided with 0.67 for x0.67 c-mount adapter
  % set these values as your setting. 
  scalefactor0 = 1;
  scalefactor100 = 9.637313/100;
  scalefactor63 = 15.30746/100;
  scalefactor40 = 18.98758/100;
% Octave and Matlab may need different code. So, here is Octave/Matlab distinction code
% This code is also product of strugle with differece between two environment. If its octave....
versionstr = version;
versionnum = str2num(versionstr(1));
% Octave2.9.9: environment == 1, Matlab6.5: environment == 0.
if(versionnum < 6)
  environment = 1;
else
  environment = 0;
end

  % suggestion for unexpected inputs.
  if(nargin < 1)
    error ('require samplename');
  elseif(nargin > 5)
    usage ('punctaanalyser(samplename, thresholdcoeff, fixthreshold, mode, objective)');
  end
  
% initial setting of thresholds and other flags
  if(nargin == 1)
    thresholdcoeff = defaultthreshold;
    fixthreshold = defaultfixthreshold;
    mode = defaultmode;
    objective = defaultobjective;
  elseif(nargin == 2)
    fixthreshold = defaultfixthreshold;
    mode = defaultmode;
    objective = defaultobjective;
  elseif(nargin == 3)
    mode = defaultmode;
    objective = defaultobjective;
  elseif(nargin == 4)
    if(mode == 0)
    mode = defaultmode;
    objective = defaultobjective;
    %Matlab dont use !=. instead ~=  
    elseif(mode ~= 0)
      mode = 1;
    end
  elseif(nargin == 5)
    if(objective == 0)
      objective = 0;
    elseif(objective == 100)
      objective = 100;
    elseif(objective == 63)
      objective = 63;
    elseif(objective == 40)
      objective = 40;
    else
      error('unknown objective: enter 0, 100, 63 or 40 at fifth value');
    end
    if(mode == 0)
      mode = defaultmode;
    elseif(mode ~= 0)
      mode = 1;
    end
  end
  if(objective == 0)
      scalefactor = scalefactor0;
  elseif(objective == 100)
      scalefactor = scalefactor100;
  elseif(objective == 63)
      scalefactor = scalefactor63;
  elseif(objective == 40)
      scalefactor = scalefactor40;
  end
  samplename
  thresholdcoeff
  fixthreshold
  mode
  objective
  scalefactor

% oneddata is sum of pic data provided by preprocessor.
% write the name of data that you want to analyse.
% filewoextension = ('wt-1')
% When you want to process all data with batch script, comment out above and use bellow instead.
%filewoextension = (imagename)
filename = strcat(samplename, '.txt');
% windows octave2.9/3.0 has bag. Cannot read tab delimited file. need to find different way.
%subpic = dlmread(filename);

%%%%%%% use below instead of dlmread
% a test code to circumvent the bug that windows octave dlmread cannnot read tab delimited imagefile
fid = fopen(filename, 'r');
%read file as vector
[picvec, num] = fscanf(fid, '%d');
fclose(fid);
fid = fopen(filename, 'rb');
% read as binary data.
binarypic = fread(fid);
fclose(fid);

% find return code 10
nr = length(find(binarypic==10));

width = num/nr;

subpic = reshape(picvec, width, nr)';
%%%%%%%

% Dont depict images if batch mode = 1
if(mode == 0)
%imagedata = imread(strcat(samplename, '.jpg'));
%figure(1);
%image(imagedata)
figure(10);
% mat2gray doesnt work in matlab, and image(subpic) gives ugry and (heavy response?) on octave/Mac.
if(environment == 1)
imshow(mat2gray(subpic))
else
image(subpic)
end
end

oneDdata = sum(subpic);
[dummy, imagewidth] = size (oneDdata);


%%%%%%%%%% making matrix for initial derivative process.
% 1st row; raw data
rawmatrix = oneDdata;
% 2nd row; derivative. To fill the last column, add 0 at the end.
rawmatrix(2,:) = [diff(rawmatrix(1,:)),0];
% 3rd row; derivative of reverse direction.To fill the 1st column, add 0 at the head.
rawmatrix(3,:) = [0, diff(rawmatrix(1,:))*-1];
% 4th row; multiplying 2nd*3rd.points that >=0 are critical points.(peak or bottm).
rawmatrix(4,:) = rawmatrix(2,:).*rawmatrix(3,:);
% 5th row;  addition 2nd+3rd. points that < 0 are peak, >0 are bottom.
%  = 0 have two possibility. one is start and end of picture. another one is flat place.
rawmatrix(5,:) = rawmatrix(2,:)+rawmatrix(3,:);

%%%%%%%%% making matrix of critical points.
% 1st row; x location
criticalmatrix = find(rawmatrix(4,:) >= 0);
% 2nd row; pixel value (bright ness)
criticalmatrix(2,:) = rawmatrix(1,criticalmatrix(1,:));
% 3rd row; the 5th row of rawmatrix.minus or plus sign indicating peak or bottm flat
criticalmatrix(3,:) = rawmatrix(5,criticalmatrix(1,:));
% 4th row; derivative of 2nd row. This is used to cut off under the threthold. see 6th row
criticalmatrix(4,:) = [diff(criticalmatrix(2,:)), 0];
% 5th row; derivative of reverse direction.
criticalmatrix(5,:) = [0, diff(criticalmatrix(2,:))*-1];
% 6th row; larger abs of 4th or 5th row. this is used to define pseudo (noise) or true critical points.
criticalmatrix(6,:) =max([abs(criticalmatrix(4,:)); abs(criticalmatrix(5,:))]);

%% Here is threshold setting.  mean oneDdata are used as standard.
% This may not good definition. fixed value would be better?
%threshold =mean(oneDdata)/100*20;
threshold =mean(oneDdata)/100*thresholdcoeff;

%%%%%%%% matrix of ture critical points
% 1st row; index of the criticalmatrix where the 6th row (absolutevalue) is larger than thresholod
truematrix = find(criticalmatrix(6,:) > threshold);
% 2nd row; location of rawdata
%size(oneDdata)
%size(truematrix)
%size(criticalmatrix)

truematrix(2,:) = criticalmatrix(1,truematrix(1,:));
% 3rd row; pixel value (bright ness)
truematrix(3,:) = criticalmatrix(2,truematrix(1,:));
% 4th row; the 3rd row of criticalmatrix.minus or plus sign indicating peak or bottm or flat point.
truematrix(4,:) = criticalmatrix(3,truematrix(1,:));

%%%%%%%% eliminating truecriticalpoints if its at start or end picture.
% for distance and linear density, use peakmatrix instead intervalmatrix
if (truematrix(2,1) == 1)
truematrix(:,1) = [];
end
[dummy, tmxcolumns] = size(truematrix);
if (truematrix(2,tmxcolumns) == imagewidth)
truematrix(:,tmxcolumns) = [];
end
[dummy, tmxcolumns] = size(truematrix);

%%%%%%%% bottom and peak matrix which has x-location
bottom = find(truematrix(4,:) >= 0);
[dummy, bottomlength] = size(bottom);
bottommatrix = truematrix(2,bottom);
peak = find(truematrix(4,:) < 0);
[dummy, peaklength] = size(peak);
peakmatrix = truematrix(2,peak);

% Dont depict images if batch mode = 1
if(mode == 0)
figure(30);
title('horizontal view')
plot(oneDdata, 'r')
hold on
%plot(criticalmatrix(1,:), criticalmatrix(2,:),'gx')
plot(bottommatrix, truematrix(3,bottom) , 'gx')
plot(peakmatrix, truematrix(3,peak) , 'b+')
hold off
end

%%%%%%% matrix of interval defined as between two bottoms
% 1st row is the start of the interval
intervalmatrix = bottommatrix(1:bottomlength-1);
% 2nd row is the end of the interval
intervalmatrix(2,:) = bottommatrix(2:bottomlength);
% 3rd row is the length of the interval
intervalmatrix(3,:) = intervalmatrix(2,:) - intervalmatrix(1,:);
% prepare 4th row for strage how many truepeaks are included within the interval at next step.
intervalmatrix(4,:) = 0;
[dummy, intervallength] = size(intervalmatrix);
for a=1:intervallength
[dummy, intervalmatrix(4,a)] = size(find((peakmatrix > intervalmatrix(1,a)) & (peakmatrix<intervalmatrix(2,a))));
end

% This data is pending. old version has it, but seems useless.
% 5th row; indicating if the interval contat with '0' gap.
% x=1:intervallength-1;
% intervalmatrix(5,x) = intervalmatrix(4,x) .* intervalmatrix(4,x+1);
% If this value >0, it means that both of the interval and the next one have peak(s) and the gap must be zero.

% 5th row has sum of intensity. punctavolume.
for a=1:intervallength
intervalmatrix(5,a) =  sum(oneDdata(intervalmatrix(1,a):intervalmatrix(2,a))) - min(oneDdata)*intervalmatrix(3,a);
% fixed 091026
%intervalmatrix(5,a) =  sum(oneDdata(intervalmatrix(1,a):intervalmatrix(2,a)) - min(oneDdata)*intervalmatrix(3,a));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing output data
% number of puncta
pcount = peaklength
% pwidths: widths of interval which have peak(s)
pwidths = intervalmatrix(3, find(intervalmatrix(4,:)>0)) .* scalefactor;
meanpwidths = mean(pwidths)
% pgaps: width of interval which dont have peak.
pgaps = intervalmatrix(3, find(intervalmatrix(4,:)==0))  .* scalefactor;
meanpgaps = mean(pgaps)
% pdist: distance from peak to peak
x=1:peaklength-1;
pdist = (peakmatrix(x+1)-peakmatrix(x)) .* scalefactor;
meanpdist = mean(pdist)
% clusterpeaks : How many peaks are included in cluster
clusterpeaks = intervalmatrix(4,find(intervalmatrix(4,:)>=2));
% intensity of peaks
pintens = truematrix(3,peak);
medianintens = median(pintens);
% normalise peak intensity with median peak intensity
% This is not used actuary.
normintens = pintens / median(truematrix(3,bottom)) * 100;
% pvolumes: volume of puncta.
pvolumes = intervalmatrix(5, find(intervalmatrix(4,:)>0));
meanpvolumes = mean(pvolumes)
% plineardensity: puncta linear density
plineardensity = pcount/imagewidth / scalefactor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% another mothod for puncta width etc.
% pix above fixthreshold as signal to measure the width.
%%This definition may not reflect actual puncta. But practically it might be work. just try.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixthreshold = 200;
% pucta must start after xmin and end before xmax. These value will bed used as end processing.
xmin = min(find(oneDdata < fixthreshold));
xmax = max(find(oneDdata < fixthreshold));

% Always oneDdatatrimmed start with less than fixthreshold and end less than fixthreshold.
oneDdatatrimmed = oneDdata(xmin:xmax);

% To detect start and of punctawidth, index which included in both above threshold and,
% less threshold moved plus and minus 1 are found.
% The starts points are storaged in the 1st row of widthmodel
% The end points are storaged in the 2nd row of widthmodel
widthmodel = [];
widthmodel = intersect (find(oneDdatatrimmed >= fixthreshold) , find(oneDdatatrimmed < fixthreshold) + 1);
widthmodel(2,:) = intersect (find(oneDdatatrimmed >= fixthreshold) , find(oneDdatatrimmed < fixthreshold) - 1);
% 3rd row: length, width what you want to know
widthmodel(3,:) = widthmodel(2,:) - widthmodel(1,:) + 1;
% 4th row: volume, accumuration of intensity
[row, column] = size(widthmodel);
for a=1:column
widthmodel(4,a) = sum(oneDdatatrimmed(widthmodel(1,a):widthmodel(2,a)));
end

fixwidth = widthmodel(3,:)  .* scalefactor;
meanfixwidth = mean(fixwidth)
fixvol = widthmodel(4,:);
meanfixvol = mean(fixvol)

% Dont depict images if batch mode = 1
if(mode == 0)
figure(40);
title('widthtes')
hold off
plot(oneDdatatrimmed, 'r')
hold on
plot(widthmodel(1,:), oneDdatatrimmed(widthmodel(1,:)), 'g+')
plot(widthmodel(2,:), oneDdatatrimmed(widthmodel(2,:)), 'kx')
%while a = 1:column
%plot((widthmodel(1,a):widthmodel(2,a)), oneDdatatrimmed(widthmodel(1,a):widthmodel(2,a)), 'b^');
%end
hold off
end

%%%%%% fixgap method silmilar with above
% gap must start after xmin and end before xmax. These value will be used as end processing.
xmin = min(find(oneDdata >= fixthreshold));
xmax = max(find(oneDdata >= fixthreshold));

% Always oneDdatatrimmed start with less than fixthreshold and end less than fixthreshold.
oneDdatatrimmed = oneDdata(xmin:xmax);

% To detect start and end of fixgap, index which included in both above threshold and,
% less threshold moved minus and plus 1 are found.
% The starts points are storaged in the 1st row of gapmodel
% The end points are storaged in the 2nd row of gapmodel
gapmodel = [];
gapmodel = intersect (find(oneDdatatrimmed >= fixthreshold) , find(oneDdatatrimmed < fixthreshold) - 1);
gapmodel(2,:) = intersect (find(oneDdatatrimmed >= fixthreshold) , find(oneDdatatrimmed < fixthreshold) + 1);
% 3rd row: length, width what you want to know
gapmodel(3,:) = gapmodel(2,:) - gapmodel(1,:) ;

fixgap = gapmodel(3,:)  .* scalefactor;
meanfixgap = mean(fixgap)

[dummy, fixwidthcount] = size(fixwidth);
[dummy, fixgapcount] = size(fixgap);
% fixwidthlineardensity: fixwidth linear density
fixwidthlineardensity = fixwidthcount/imagewidth / scalefactor
% fixgaplineardensity: fixgap linear density
fixgaplineardensity = fixgapcount/imagewidth / scalefactor


end
