clear all

MasterSpreadSheetID = '1D_4OouMg8elj_djLPVKltK7Yq_Urgayx4KI2zEomAEE';
diffusionfile = 'C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Observed_Data\Diffusion_Coefficient_Data.csv';


[mastersheet, edaddress, pathname, pathnameA, pathnameB, pathnameC, pathnameF, pathnameG, filenameD, Final_Results, total] = pulldata(MasterSpreadSheetID);
[difflist] = dcoeffpull(mastersheet);
trunctime = [];
cutoff = .72; %cut off solvent height growth, cm

set(groot,'defaultLineLineWidth',2.0)

%% Creates an array of the experimental data, and inputs initial conditions

for exp = 1:46
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

%%

% Truncating the solvent readings

solventd = strunc(solventd,cutoff);


% Rounding to the timestep
dt = 20; %timestep in seconds, standard is 10 seconds
for i = 1:length(solventd)
    solventd(i,2) = dt * round(solventd(i,2)/dt);
end

%% Initial Conditions, adjust as needed to fit the observed data

%Input the initial setting of the run
diffusionuncertainty = 2; %diffusion uncertainty. 0 is exp, 1 is low, 2 is high
calibrationrun = 0; %is this a calibration run?


columnheight = solventd(1,1); %height of column in cm
radius = .35; %.3313; %radius of column in cm
m0 = 0; %inital mols of antisolvent in bins
bins = 10; %number of bins
tduration = max(solventd(:,2)); %length of experiment in seconds
d = difflist(exp,:) * 10000; % given diffusion constant in m^2/s, converted to cm^2/s
dorig = d; %Keeps track of original diffusion constants
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


%% Crystallization Pull

[xtaltime,xtalheight,xtaldata] = crystallizationpull(fribble,calibrationrun,rxnid);

%rounds xtaltime to correct dt
if xtaldata == "Yes"
    xtaltime = dt * round(xtaltime/dt);
    xtalindex = xtaltime/dt;
end
    
%%  Starting Conditions
[antimols,orgmols,inorgmols,solvmols,FAHmols,solventvolume,orgvol,inorgvol,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind] = StartingCalculationsV3(si,fahi,solvent,orgid,orgsol,inorgsol);

parameters = 'C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Tmp\parameters.mat';
save(parameters, 'rxnid','solventd','columnheight','radius','m0','bins','dt','tduration','d','mva','cs0','co0','antimols','solvmols','inorgmols','orgmols','FAHmols','htov','solventord','organicord','solventconv','organicconv','antiind','solvind','orgind','FAHind','inorgind','vtotal')

%% Creates Folder for Experiment to Output Data, And creates Final Path

%filenameA = fullfile(pathnameA, rxnid + '_Best_Fits.png');
%filenameB = fullfile(pathnameB, rxnid + '_Concentration_Profile.png');
%filenameC = fullfile(pathnameC, rxnid + '_Crystallization_Profile.png');

%% Visualization Loop Start, 3 loop, one for each diffusion coefficient, edits Final Path accordingly

for visloop = 1:3
    
% Final Fit Plot


%Adjusting diffusion for uncertainty, 1 is expected, 2 is low, 3 is high, edits Final Path accordingly

if visloop == 1
    d = dorig(1,1);
    disp("diffusion expected run")
    filenameA = fullfile(pathnameA, rxnid + ' Expected_Best_Fits.png');
    filenameB = fullfile(pathnameB, rxnid + ' Expected_Concentration_Profile.png');
    filenameC = fullfile(pathnameC, rxnid + ' Crystallization_Profile.png');
    filenameE = fullfile(pathname, 'MinimumsExp\' + rxnid + '.mat');
    filenameF = fullfile(pathnameF, rxnid + ' Final_Fits.png');
    filenameG = fullfile(pathnameG, rxnid + ' ExpectedColumnVisualization.png');
    bound = "Expected Diffusion Coefficient";
elseif visloop == 2
    d = dorig(1,1)- dorig(1,2);
    if d < 0
        d = 0;
    end
    disp("diffusion low uncertainty run")
    filenameA = fullfile(pathnameA, rxnid + ' LowerBound_Best_Fits.png');
    filenameB = fullfile(pathnameB, rxnid + ' LowerBound_Concentration_Profile.png');
    filenameC = fullfile(pathnameC, rxnid + ' Crystallization_Profile.png');
    filenameE = fullfile(pathname, 'MinimumsLow\' + rxnid + '.mat');
    filenameF = fullfile(pathnameF, rxnid + ' Final_Fits.png');
    filenameG = fullfile(pathnameG, rxnid + ' LowerBoundColumnVisualization.png');
    bound = "Lower Bound Diffusion Coefficient";
else
    d = dorig(1,1) + dorig(1,2);
    disp("diffusion high uncertainty run")
    filenameA = fullfile(pathnameA, rxnid + ' UpperBound_Best_Fits.png');
    filenameB = fullfile(pathnameB, rxnid + ' UpperBound_Concentration_Profile.png');
    filenameC = fullfile(pathnameC, rxnid + ' Crystallization_Profile.png');
    filenameE = fullfile(pathname, 'MinimumsHigh\' + rxnid + '.mat');
    filenameF = fullfile(pathnameF, rxnid + ' Final_Fits.png');
    filenameG = fullfile(pathnameG, rxnid + ' UpperBoundColumnVisualization.png');
    bound = "Upper Bound Diffusion Coefficient";
end



load(filenameE)
ccat = [];
timecheck = randi(tduration/dt);
basepos = .01; %Some random position to feed the model query
minforpc = 4;
sminforpc = num2str(minforpc);

for pling = 1:length(mins)
[t,hprofile,cprofile,columnheightout,~,~,bincentout] = DiffusionAllSpeciesV2Visualization(columnheight,radius,m0,bins,dt,tduration,d,mva,cs0,co0,mins(pling,1),mins(pling,2),antimols,solvmols,inorgmols,orgmols,FAHmols,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind,vtotal,basepos,timecheck,0);
% Making an array of the observed timepoints in the model data
modeldata = cat(2,transpose(columnheightout(1:end-1)),transpose(t));
d2comp = zeros(length(solventd),1);


for i = 1:length(solventd)
    d2comp(i) = modeldata(find(modeldata(:,2)==solventd(i,2)),1);
end
d2comp = cat(2,d2comp,solventd);

%Making an array of all of the concentration profiles
ccat = cat(1,ccat,cprofile);

skin = num2str(mins(pling,1));
skout = num2str(mins(pling,2));
ssse = num2str(mins(pling,3));
fitnum = num2str(pling);
stimecheck = num2str(timecheck);
% Plotting the fit for the 10 minimums
fig1 = figure(1);
subplot(5,2,pling)
%sgtitle(['Model and Experimental Solution Height over Time for ' char(rxnid) bound],'FontWeight','bold')
plot(d2comp(:,3),d2comp(:,1))
hold on
plot(d2comp(:,3),d2comp(:,2))
legend({['Model SSE: ' ssse],'Experimental'},'Location','southeast')
title(['Fit ' fitnum] )
xlabel('Seconds','FontWeight','bold')
ylabel('Height (cm)','FontWeight','bold')
hold off
ax = gca;
ax.LineWidth = 2;
set(gcf,'Position',get(0,'Screensize'));
%saveas(gca,filenameA) %saves figure


%Takes Best Fit
if pling == 1
    BestFitsPlots = figure(5);
    sgtitle(['Best Fits for ' char(rxnid)],'FontWeight','bold')
    subplot(1,3,visloop)
    plot(d2comp(:,3),d2comp(:,1))
    hold on
    plot(d2comp(:,3),d2comp(:,2))
    legend({['Model | Kin: ' skin ' mols/sec | Kout: ' skout ' mols/sec | SSE: ' ssse],'Experimental'},'Location','southeast')
    if visloop == 1
        title('Expected Diffusion Constant','FontWeight','bold')
    elseif visloop == 2
        title('Lower Bound Diffusion Constant','FontWeight','bold')
    else
        title('Upper Bound Diffusion Constant','FontWeight','bold')
    end
    xlabel('Seconds','FontWeight','bold')
    ylabel('Height (cm)','FontWeight','bold')
    hold off
    ax = gca;
    ax.LineWidth = 2;
    set(gcf,'Position',get(0,'Screensize'));
    %saveas(gca,filenameB) %saves figure
end

% Saves D2COMP FOR LATER PLOTTING
if pling == 1 && visloop == 1
    d2comptobeplotted = d2comp;
end

%Random Concentration Profile Check
if pling <= minforpc
    fig2 = figure(2);
    hold on
    title(['Fits 1-' sminforpc ' Superimposed Concentration Profile Check at Random Time, ' stimecheck ' seconds'])
    plot(bincentout(:,timecheck),cprofile(:,timecheck))
    xlabel('Height (cm)','FontWeight','bold')
    ylabel('Mols/L','FontWeight','bold')
    hold off
    ax = gca;
    ax.LineWidth = 2;
    set(gcf,'Position',get(0,'Screensize'));
    %saveas(gca,filenameB) %saves figure
end




end

%Final Edits and saves the figures 
%Fit figure
ax = gca;
ax.LineWidth = 2;
set(gcf,'Position',get(0,'Screensize'));
supersizeme(1.35) %Adam Danz (2021). supersizeme (https://www.mathworks.com/matlabcentral/fileexchange/67644-supersizeme), MATLAB Central File Exchange. Retrieved February 11, 2021.
saveas(gca,filenameA) %saves figure
close(fig1)
%Profile figure
ax = gca;
ax.LineWidth = 2;
set(gcf,'Position',get(0,'Screensize'));
supersizeme(1.5) %Adam Danz (2021). supersizeme (https://www.mathworks.com/matlabcentral/fileexchange/67644-supersizeme), MATLAB Central File Exchange. Retrieved February 11, 2021.
saveas(gca,filenameB) %saves figure
close(fig2)

%close all

%% Final Visualization and Results
if xtaldata == "Yes"
[t,hprofile,cprofile,columnheightout,~,~,bincentout,dcmout,fahout,orgout,inorgout,solvout] = DiffusionAllSpeciesV2Visualization(columnheight,radius,m0,bins,dt,tduration,d,mva,cs0,co0,mins(1,1),mins(1,2),antimols,solvmols,inorgmols,orgmols,FAHmols,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind,vtotal,xtalheight,xtalindex,0);


% Making an array of the observed timepoints in the model data
modeldata = cat(2,transpose(columnheightout(1:end-1)),transpose(t));
d2comp = zeros(length(solventd),1);

for i = 1:length(solventd)
    d2comp(i) = modeldata(find(modeldata(:,2)==solventd(i,2)),1);
end
d2comp = cat(2,d2comp,solventd);

%Records the data for each run given the diffusion constant and minimums
if visloop == 1
    expdcm = dcmout;
    expfah = fahout;
    exporg = orgout;
    expinorg = inorgout;
    expsolv = solvout;
    expheights = bincentout(:,xtalindex);
    exptotheight = sum(hprofile(:,xtalindex));
    expkin = mins(1,1);
    expkout = mins(1,2);
    expseedkin = mins (1,4);
    expseedkout = min(1,5);
    expsse = mins(1,3);
    expdiffconstant = d; %cm^2/sec
    expdcmint = LinearInterpV2(expheights,expdcm,xtalheight);
    expfahint = LinearInterpV2(expheights,expfah,xtalheight);
    exporgint = LinearInterpV2(expheights,exporg,xtalheight);
    expinorgint = LinearInterpV2(expheights,expinorg,xtalheight);
    expsolvint = LinearInterpV2(expheights,expsolv,xtalheight);
elseif visloop == 2
    lowdcm = dcmout;
    lowfah = fahout;
    loworg = orgout;
    lowinorg = inorgout;
    lowsolv = solvout;
    lowheights = bincentout(:,xtalindex);
    lowtotheight = sum(hprofile(:,xtalindex));
    lowkin = mins(1,1);
    lowkout = mins(1,2);
    lowseedkin = mins (1,4);
    lowseedkout = min(1,5);
    lowsse = mins(1,3);
    lowdiffconstant = d; %cm^2/sec
    lowdcmint = LinearInterpV2(lowheights,lowdcm,xtalheight);
    lowfahint = LinearInterpV2(lowheights,lowfah,xtalheight);
    loworgint = LinearInterpV2(lowheights,loworg,xtalheight);
    lowinorgint = LinearInterpV2(lowheights,lowinorg,xtalheight);
    lowsolvint = LinearInterpV2(lowheights,lowsolv,xtalheight);
else
    highdcm = dcmout;
    highfah = fahout;
    highorg = orgout;
    highinorg = inorgout;
    highsolv = solvout;
    highheights = bincentout(:,xtalindex);
    hightotheight = sum(hprofile(:,xtalindex));
    highkin = mins(1,1);
    highkout = mins(1,2);
    highseedkin = mins (1,4);
    highseedkout = min(1,5);
    highsse = mins(1,3);
    highdiffconstant = d; %cm^2/sec
    highdcmint = LinearInterpV2(highheights,highdcm,xtalheight);
    highfahint = LinearInterpV2(highheights,highfah,xtalheight);
    highorgint = LinearInterpV2(highheights,highorg,xtalheight);
    highinorgint = LinearInterpV2(highheights,highinorg,xtalheight);
    highsolvint = LinearInterpV2(highheights,highsolv,xtalheight);
end

else %If there is no crystallization, we need to fill out values for the table
xtalheight = solventd(length(solventd)-1,1);%top %solventd(length(solventd),1)/2; halfway
xtaltime = solventd(length(solventd),1);
xtalindex = solventd(length(solventd),2)/dt;
[t,hprofile,cprofile,columnheightout,~,~,bincentout,dcmout,fahout,orgout,inorgout,solvout] = DiffusionAllSpeciesV2Visualization(columnheight,radius,m0,bins,dt,tduration,d,mva,cs0,co0,mins(1,1),mins(1,2),antimols,solvmols,inorgmols,orgmols,FAHmols,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind,vtotal,xtalheight,xtalindex,0);

if visloop == 1
    expdcm = dcmout;
    expfah = fahout;
    exporg = orgout;
    expinorg = inorgout;
    expsolv = solvout;
    expheights = bincentout(:,xtalindex);
    exptotheight = sum(hprofile(:,xtalindex));
    expkin = mins(1,1);
    expkout = mins(1,2);
    expseedkin = mins (1,4);
    expseedkout = min(1,5);
    expsse = mins(1,3);
    expdiffconstant = d; %cm^2/sec
    expdcmint = LinearInterpV2(expheights,expdcm,xtalheight);
    expfahint = LinearInterpV2(expheights,expfah,xtalheight);
    exporgint = LinearInterpV2(expheights,exporg,xtalheight);
    expinorgint = LinearInterpV2(expheights,expinorg,xtalheight);
    expsolvint = LinearInterpV2(expheights,expsolv,xtalheight);
elseif visloop == 2
    lowdcm = dcmout;
    lowfah = fahout;
    loworg = orgout;
    lowinorg = inorgout;
    lowsolv = solvout;
    lowheights = bincentout(:,xtalindex);
    lowtotheight = sum(hprofile(:,xtalindex));
    lowkin = mins(1,1);
    lowkout = mins(1,2);
    lowseedkin = mins (1,4);
    lowseedkout = min(1,5);
    lowsse = mins(1,3);
    lowdiffconstant = d; %cm^2/sec
    lowdcmint = LinearInterpV2(lowheights,lowdcm,xtalheight);
    lowfahint = LinearInterpV2(lowheights,lowfah,xtalheight);
    loworgint = LinearInterpV2(lowheights,loworg,xtalheight);
    lowinorgint = LinearInterpV2(lowheights,lowinorg,xtalheight);
    lowsolvint = LinearInterpV2(lowheights,lowsolv,xtalheight);
else
    highdcm = dcmout;
    highfah = fahout;
    highorg = orgout;
    highinorg = inorgout;
    highsolv = solvout;
    highheights = bincentout(:,xtalindex);
    hightotheight = sum(hprofile(:,xtalindex));
    highkin = mins(1,1);
    highkout = mins(1,2);
    highseedkin = mins (1,4);
    highseedkout = min(1,5);
    highsse = mins(1,3);
    highdiffconstant = d; %cm^2/sec
    highdcmint = LinearInterpV2(highheights,highdcm,xtalheight);
    highfahint = LinearInterpV2(highheights,highfah,xtalheight);
    highorgint = LinearInterpV2(highheights,highorg,xtalheight);
    highinorgint = LinearInterpV2(highheights,highinorg,xtalheight);
    highsolvint = LinearInterpV2(highheights,highsolv,xtalheight);
end

end

%Gets the column
[t,hprofile,cprofile,columnheightout,~,~,bincentout,dcmout,fahout,orgout,inorgout,solvout,visualization] = DiffusionAllSpeciesV2Visualization(columnheight,radius,m0,bins,dt,tduration,d,mva,cs0,co0,mins(1,1),mins(1,2),antimols,solvmols,inorgmols,orgmols,FAHmols,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind,vtotal,xtalheight,xtalindex,1);
set(gcf,'Position',get(0,'Screensize'));
supersizeme(1.5) %Adam Danz (2021). supersizeme (https://www.mathworks.com/matlabcentral/fileexchange/67644-supersizeme), MATLAB Central File Exchange. Retrieved February 11, 2021.
saveas(gca,filenameG)
close

end % end of visloop, data for all three diffusion conditions has been collected at this point


% Formats the uncertainty window

if xtaldata == "Yes"
    %Uncertainty Range
    %DCM
    dcmrange = cellstr(strcat('[',num2str(max([expdcmint,lowdcmint,highdcmint])),',',num2str(min([expdcmint,lowdcmint,highdcmint])),']'));
    %FAH
    fahrange = cellstr(strcat('[',num2str(max([expfahint,lowfahint,highfahint])),',',num2str(min([expfahint,lowfahint,highfahint])),']'));
    %ORG
    orgrange = cellstr(strcat('[',num2str(max([exporgint,loworgint,highorgint])),',',num2str(min([exporgint,loworgint,highorgint])),']'));
    %INORG
    inorgrange = cellstr(strcat('[',num2str(max([expinorgint,lowinorgint,highinorgint])),',',num2str(min([expinorgint,lowinorgint,highinorgint])),']'));
    %SOLV
    solvrange = cellstr(strcat('[',num2str(max([expsolvint,lowsolvint,highsolvint])),',',num2str(min([expsolvint,lowsolvint,highsolvint])),']'));
else
    %DCM
    dcmrange = cellstr(strcat('[',num2str(max([expdcmint,lowdcmint,highdcmint])),',',num2str(min([expdcmint,lowdcmint,highdcmint])),']'));
    %FAH
    fahrange = cellstr(strcat('[',num2str(max([expfahint,lowfahint,highfahint])),',',num2str(min([expfahint,lowfahint,highfahint])),']'));
    %ORG
    orgrange = cellstr(strcat('[',num2str(max([exporgint,loworgint,highorgint])),',',num2str(min([exporgint,loworgint,highorgint])),']'));
    %INORG
    inorgrange = cellstr(strcat('[',num2str(max([expinorgint,lowinorgint,highinorgint])),',',num2str(min([expinorgint,lowinorgint,highinorgint])),']'));
    %SOLV
    solvrange = cellstr(strcat('[',num2str(max([expsolvint,lowsolvint,highsolvint])),',',num2str(min([expsolvint,lowsolvint,highsolvint])),']'));
end

%Plotting Concentration Profile for All Species if Crystallization


concplot = figure(3);
if xtaldata == "Yes"
    sgtitle(['Concentration Profiles for ' char(rxnid) ' at Crystallization'],'FontWeight','bold')
else
    sgtitle(['No Crystallization for ' char(rxnid) ', Concentration Profile at End of Expt, Top of Column'],'FontWeight','bold')
end

subplot(1,2,1)
plot(d2comptobeplotted(:,3),d2comptobeplotted(:,1))
hold on
plot(d2comptobeplotted(:,3),d2comptobeplotted(:,2))
legend({['Model | Kin: ' skin ' mols/sec | Kout: ' skout ' mols/sec | SSE: ' ssse],'Experimental'},'FontWeight','bold','Location','southeast')
title('Expected Diffusion Constant Model to Experimental Data','FontWeight','bold')
xlabel('Seconds','FontWeight','bold')
ylabel('Height (cm)','FontWeight','bold')
hold off
ax = gca;
ax.LineWidth = 2;
set(gcf,'Position',get(0,'Screensize'));
    
subplot(1,2,2)
hold on
plot(expheights,expdcm)
plot(lowheights,lowdcm)
plot(highheights,highdcm)
axis([0 min([hightotheight,lowtotheight,exptotheight]) 0 inf])
xlabel("Height (cm)",'FontWeight','bold')
ylabel("Concentration (mols/liter)",'FontWeight','bold')
xline(xtalheight,'-.','Crystallization','FontWeight','bold','LineWidth',4)
legend({['Expected Diffusion Coefficient ' num2str(expdiffconstant) ' cm^2/sec'],['Lower Bound Diffusion Coefficient ' num2str(lowdiffconstant) ' cm^2/sec'],['Upper Bound Diffusion Coefficient ' num2str(highdiffconstant) ' cm^2/sec'],['M=' char(dcmrange) ' mols/L' newline 'Expected M=' num2str(expdcmint) ' mols/L']},'FontWeight','bold','Location','southwest')
ax = gca;
ax.LineWidth = 2;
hold off
title("DCM Concentration (mols/liter) per Height (cm)")




%Crystallization Figure
ax = gca;
ax.LineWidth = 3;
set(gcf,'Position',get(0,'Screensize'));
supersizeme(1.7) %Adam Danz (2021). supersizeme (https://www.mathworks.com/matlabcentral/fileexchange/67644-supersizeme), MATLAB Central File Exchange. Retrieved February 11, 2021.
saveas(gca,filenameC)
close

%Final Fits Figures
ax = gca;
ax.LineWidth = 2;
set(gcf,'Position',get(0,'Screensize'));
supersizeme(1.75) %Adam Danz (2021). supersizeme (https://www.mathworks.com/matlabcentral/fileexchange/67644-supersizeme), MATLAB Central File Exchange. Retrieved February 11, 2021.
saveas(gca,filenameF)

Final_Results = vertcat(Final_Results, table(rxnid,solvent,amine,expdiffconstant,expkin,expkout,expsse,xtaldata,xtaltime,xtalheight,expdcmint,dcmrange,expsolvint,solvrange,expfahint,fahrange,expinorgint,inorgrange,exporgint,orgrange,'VariableNames',{'Reaction ID','Solvent','Amine','Expected Diffusion Coefficient','Best Kin, Expected Diffusion Coeff, Mols/sec','Best Kout, Expected Diffusion Coeff, Mols/sec','Best Fit SSE, Expected Diffusion Coeff','Crystallization Data?','Crystallization Time, Seconds','Crystallization Height, cm','Expected DCM Concentration, Mols/L','DCM Concentration Uncertainty Range','Expected Solvent Concentration, Mols/L','Solvent Concentration Uncertainty Range','Expected FAH Concentration, Mols/L','FAH Concentration Uncertainty Range','Expected Inorganic Concentration, Mols/L','Inorganic Concentration Uncertainty Range','Expected Organic Concentration, Mols/L','Organic Concentration Uncertainty Range'}));

close all

end %end of total loop

writetable(Final_Results, filenameD)

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
function [mastersheet, edaddress, pathname, pathnameA, pathnameB, pathnameC, pathnameF, pathnameG, filenameD, Final_Results, total] = pulldata(MasterSpreadSheetID)
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
    pathnameF = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\Final_Fits_Plot\');
    pathnameG = fileparts('C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Results\Column_Visualization\');
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

%% Gets Crystallization Location and Time

function [xtaltime,xtalheight,xtaldata] = crystallizationpull(fribble,calibrationrun,rxnid)
wheressolv = find(contains(fribble(1,:),'measurement_name') == 1);
wheresxheight = find(contains(fribble(1,:),'height_cm_of_mean_crystal_vertex_from_bottom') == 1);
wheresxtal = find(contains(fribble(1,:),'crystal_present') == 1);
wherest = find(contains(fribble(1,:),'elapsed_minutes') == 1);

issolv = find(strcmp(fribble(:,wheressolv),'solvent') == 1);



fribble(1,9) = {''};
hcrys = strcmp(fribble(:,9),'');

if sum(str2double(fribble(2:end,wheresxtal))) == 0
    xtalad = "";
    if calibrationrun == 1
        xtalad = "calibration";
    end
else
    rowostuff = intersect(issolv, find(hcrys == 0));
    xtalad = "crystallization";
end


if xtalad == "crystallization"
    disp 'crystallization data detected'
    xtaldata = "Yes";
    xtaltime = round(str2double(fribble(rowostuff,wherest)) * 60); %% In seconds
    xtalheight = str2double(fribble(rowostuff,wheresxheight));
elseif xtalad == "calibration"
    disp 'Calibration Run, random time and height selected'
    xtaldata = "Test Check";
    xtaltime = randi(length(t));
    xtalheight = rand * columnheightout(1,xtaltime);
else
    disp 'No crystallization data for'
    disp (rxnid)
    xtaldata = "None, Crystallization Concentration Profile Taken at End of Exp, Middle of Column";
    xtaltime = "NA";
    xtalheight = "NA";
end

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

