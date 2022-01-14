clear all

MasterSpreadSheetID = '1D_4OouMg8elj_djLPVKltK7Yq_Urgayx4KI2zEomAEE';
diffusionfile = 'C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Observed_Data\Diffusion_Coefficient_Data.csv';

[mastersheet, edaddress, pathname, pathnameA, pathnameB, pathnameC, filenameD, Final_Results, total] = pulldata(MasterSpreadSheetID);
[difflist] = dcoeffpull(mastersheet);
trunctime = [];
cutoff = .72; %cut off solvent height growth, cm


%% Creates an array of the experimental data, and inputs initial conditions

for exp = 1:1
test = string(edaddress(exp));

fribble = GetGoogleSpreadsheet(convertStringsToChars(test));
frobble = strfind(fribble,"solvent");
frabble = cell2mat(frobble);
frabble = cat(1,999,frabble);
inds = find(frabble == 1);
solventd = str2double(fribble(inds,12:2:end));
rxnid = string(mastersheet(exp+1,1));
rxnid = strrep(rxnid,'_',' ');
rxnid = strrep(rxnid,':',' ');
solvent = string(mastersheet(exp+1,find(strcmp(mastersheet(1,:),'Solvent') == 1)));
amine = string(mastersheet(exp+1,find(strcmp(mastersheet(1,:),'Amine') == 1))); % 'none' %set to none if calibration run
orgid = getoldstring(amine);
% Convert Solvent Readings to seconds from minutes
solventd(:,2) = solventd(:,2) * 60;

%Input the Results Filename ,filenameE is made in the diffusion uncertainty
%portion, ~line 92

filenameA = fullfile(pathnameA, rxnid + '_Best_Fits.png');
filenameB = fullfile(pathnameB, rxnid + '_Concentration_Profile.png');
filenameC = fullfile(pathnameC, rxnid + '_Crystallization_Profile.png');
%filenameE = fullfile(pathname, 'MinimumsExp\' + rxnid + '.mat'); 

% Truncating the solvent readings

solventd = strunc(solventd,cutoff);


% Rounding to the timestep
dt = 20; %timestep in seconds, standard is 10 seconds
for i = 1:length(solventd)
    solventd(i,2) = dt * round(solventd(i,2)/dt);
end

%% Initial Conditions, adjust as needed to fit the observed data

diffusionuncertainty = 2; %diffusion uncertainty. 0 is exp, 1 is low, 2 is high
calibrationrun = 0; %is this a calibration run?
columnheight = solventd(1,1); %height of column in cm
radius = .35; %.3313; %radius of column in cm
m0 = 0; %inital mols of antisolvent in bins
bins = 10; %number of bins
tduration = max(solventd(:,2)); %length of experiment in seconds
d = difflist(exp,:) * 10000; % given diffusion constant in m^2/s, converted to cm^2/s
mva = 64.05; %molar volume of antisolvent (DCM)
si = str2double(mastersheet(exp+1,find(strcmp(mastersheet(1,:),'Solvent vol (uL)') == 1)));
fahi = str2double(mastersheet(exp+1,find(strcmp(mastersheet(1,:),'FAH (uL)') == 1)));
invol = si + fahi;


inorgsol = str2double(mastersheet(exp+1,find(strcmp(mastersheet(1,:),'Actual In-organic [M]') == 1)));
if isempty(inorgsol) == 1
    inorgsol = 0;
end

orgsol = str2double(mastersheet(exp+1,find(strcmp(mastersheet(1,:),'Actual Organic [M]') == 1)));
if isempty(orgsol) == 1
    orgsol = 0;
end

[smolarity,fmolarity,smol,fmol,vtotal] = StartingCalculationsV2(solvent,si,fahi,invol,inorgsol,orgsol,amine); % Gives starting Concentrations of solvent, FAH

cs0 = inorgsol + orgsol;
if isempty(cs0) == 1
    cs0 = 0;
end
co0 =  smolarity + fmolarity; 



kin = .01; %condensation rate
kout = .01; %evaporation rate
k = [kin,kout];

%Adjusting diffusion for uncertainty, 1 is expected, 2 is low, 3 is high

if diffusionuncertainty == 0
    d = d(1,1);
    disp("diffusion expected run")
    filenameE = fullfile(pathname, 'MinimumsExp\' + rxnid + '.mat');
elseif diffusionuncertainty == 1
    d = d(1,1)- d(1,2);
    if d < 0
        d = 0;
    end
    disp("diffusion low uncertainty run")
    filenameE = fullfile(pathname, 'MinimumsLow\' + rxnid + '.mat');
else
    d = d(1,1) + d(1,2);
    disp("diffusion high uncertainty run")
    filenameE = fullfile(pathname, 'MinimumsHigh\' + rxnid + '.mat');
end
    

%%  Starting Conditions
[antimols,orgmols,inorgmols,solvmols,FAHmols,solventvolume,orgvol,inorgvol,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind] = StartingCalculationsV3(si,fahi,solvent,orgid,orgsol,inorgsol);

parameters = 'C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Tmp\parameters.mat';
save(parameters, 'rxnid','solventd','columnheight','radius','m0','bins','dt','tduration','d','mva','cs0','co0','antimols','solvmols','inorgmols','orgmols','FAHmols','htov','solventord','organicord','solventconv','organicconv','antiind','solvind','orgind','FAHind','inorgind','vtotal')



%% Nelder meade Fitting

completetime = tic;
results = [];
options = optimset('MaxFunEvals',100);
for i = 1:10
    if i > 5
        options = optimset('MaxFunEvals',50);
    end
k0 = [.0000000001*10^i,.0000000001*10^i]; % Initital k in and kout values. k(1) is kin, k(2) is kout
[fit,sse] = fminsearch(@DiffusionSSE3,k0,options);
results = cat(1,results,[fit,sse,k0]);
display (k0)
end
completeelapsed = toc(completetime);
%% Ordered Fits and table creation
mins = [];
for i=1:10
  [M,I] = min(results(:,3));
  mins = cat(1,mins,results(I,:));
  results(I,:) = [];
end
%%

save(filenameE,'mins')
%% Closing Figures
%close all
end
%writetable(Final_Results, filenameD)




%% FUNCTIONS

% Truncating the solvent data to the assigned cutoff
function [trssolv] = strunc(solventd,cutoff)
    sbu = [solventd(:,1) - solventd(1,1),solventd(:,2)];
    [minValue, closestIndex] = min(abs(sbu(:,1) - cutoff));
    closestValue = sbu(closestIndex);
    trssolv = solventd(1:closestIndex,:);
    %trunctime = cat(1,trunctime,[rxnid,solventd(closestIndex,:)]); %Used to create a record of all cutoff points
    
end

% Getting the info from google drive
function [mastersheet, edaddress, pathname, pathnameA, pathnameB, pathnameC, filenameD, Final_Results, total] = pulldata(MasterSpreadSheetID)
    mastersheet = GetGoogleSpreadsheet(MasterSpreadSheetID); %%Insert ID of spreadsheet
    address = mastersheet(2:end,find(strcmp(mastersheet(1,:),'Diffusion_height_3') == 1));
    address = strcat(address, '/');
    edaddress = erase(address, 'https://docs.google.com/spreadsheets/d/');
    edaddress = erase(eraseBetween(edaddress, '/','/'),'//');
    total = length(edaddress);

    % Creates empty table to store final results
    Final_Results = table();
    %Input the path of the results folder REMEMBER TO UPDATE
    pathname = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\'); % ResultsFolder
    pathnameA = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\Best_Fits\');
    pathnameB = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\Concentration_Profile\');
    pathnameC = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\Crystallization_Profile\');
    filenameD = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\Final_Results_Table\Final_Results.csv\');
end

%% Gets Diffusion Data, Calculates Average, Reformats into Array
function [difflist] = dcoeffpull(mastersheet)
    %as of 1/23/21 from KDE and lowest uncertainty, 1st element is the
    %expected, second is the uncertainty range
    dmso = [1.23601763141488E-09,2.78549420854698E-10];
    gbl = [5.26279754219915E-10,4.55652245519834E-10];
    dmf = [5.55339190057498E-10,5.65667777931015E-10];
    dmfdmso = [2.07481613378782E-10,7.50456056585554E-11];
    gbldmf = [1.08310935431657E-10,1.04227661254505E-10];
        
    exsolvlist = mastersheet(2:end,find(strcmp(mastersheet(1,:),'Solvent') == 1));
    wdmso = double(strncmpi(exsolvlist,'DMSO',10)) .* dmso;
    wgbl = double(strncmpi(exsolvlist,'GBL',10)) .* gbl;
    wdmf = double(strncmpi(exsolvlist,'DMF',10)) .* dmf;
    wdmfdmso = double(strncmpi(exsolvlist,'DMF:DMSO',10)) .* dmfdmso;
    wgbldmf = double(strncmpi(exsolvlist,'GBL:DMF',10)) .* gbldmf;

    difflist = wdmso + wgbl + wdmf + wdmfdmso + wgbldmf;
end

%% Replace Amine Name with Old String
function [oldstring] = getoldstring(newstring)

if newstring == "aep"
    oldstring = "2Pyrrolidin1ium1ylethylammoniumiodide";
elseif newstring == "chma"
    oldstring = "CyclohexylmethylammoniumIodide";
elseif newstring == "acet"
    oldstring = "AcNH3I";
elseif newstring == "ea"
    oldstring = "EtNH3I";
elseif newstring == "ma"
    oldstring = "MeNH3I";
elseif newstring == "phenea"
    oldstring = "PhEtNH3I";
elseif newstring == "1,3-dap"
    oldstring = "Propane13diammoniumIodide";
elseif newstring == "dmed"
    oldstring = "NNDimethylethane12diammoniumiodide";
elseif newstring == "dedap"
    oldstring = "NNDiethylpropane13diammoniumiodide";
elseif newstring == "dabz"
    oldstring = "Benzenediaminedihydroiodide";
end

end