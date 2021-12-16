%% Starting Chemical Properties and Array Creation
function [antimols,orgmols,inorgmols,solvmols,FAHmols,solventvolume,orgvol,inorgvol,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind] = StartingCalculationsV3(solvol,FAHvol,solventid,orgid,orgconc,inorgconc)
%Height to volume relationship (experimentally determined)
htov = .002966; % cm/microliter

%densities
gblden = 1.12; %g/ml
dmfden = .948;
dmsoden = 1.1;
gbldmfden = 1.034;
dmfdmsoden = 1.024;
fahden = 1.22;
pbi2den = 6.16;
dcmden = 1.33;

%molar masses
gblmm = 86.090; %g/mol
dmfmm = 73.095;
dmsomm = 78.13;
gbldmfmm = 79.593;
dmfdmsomm = 75.613;
fahmm = 46.03;
pbi2mm = 461.01;
dcmmm = 84.93;

%Organic Amine Density
Pden = 2.079;
Cyden = 1.505;
Acden = 2.1766;
Etden = 2.053;
Meden = 2.341;
Phden = 1.630;
Prden = 2.102;
NNeden = 2.180;
NNpden = 1.848;
Beden = 2.341;

%Organic Amine Molar Mass
Pmm = 370.02;
Cymm = 241.11;
Acmm = 185.97;
Etmm = 173.00;
Memm = 158.97;
Phmm = 249.10;
Prmm = 329.95;
NNemm = 343.98;
NNpmm = 386.06;
Bemm = 363.97;


solventarray = [gblden,dmfden,dmsoden,gbldmfden,dmfdmsoden,fahden,pbi2den,dcmden;gblmm,dmfmm,dmsomm,gbldmfmm,dmfdmsomm,fahmm,pbi2mm,dcmmm];
solventord = {'GBL','DMF','DMSO','GBL:DMF','DMF:DMSO','FAH','PBI2','DCM'};

organicarray = [Pden,Cyden,Acden,Etden,Meden,Phden,Prden,NNeden,NNpden,Beden;Pmm,Cymm,Acmm,Etmm,Memm,Phmm,Prmm,NNemm,NNpmm,Bemm];
organicord = {'2Pyrrolidin1ium1ylethylammoniumiodide','CyclohexylmethylammoniumIodide','AcNH3I','EtNH3I','MeNH3I','PhEtNH3I','Propane13diammoniumIodide','NNDimethylethane12diammoniumiodide','NNDiethylpropane13diammoniumiodide','Benzenediaminedihydroiodide'};


solventconv = (solventarray(2,:) .* (solventarray(1,:).^-1))./ 1000; %converts moles to liters
organicconv = (organicarray(2,:) .* (organicarray(1,:).^-1))./ 1000; %converts moles to liters

%% Calculations
%anticonc = anticonc; % [M]
%solvconc = solvconc; % [M]
orgconc = orgconc; % [M]
%FAHconc = FAHconc; % [M]
inorgconc = inorgconc; % [M]

antiid = 'DCM';
solventid = solventid;
orgid = orgid;
FAHid = 'FAH';
inorgid = 'PBI2';

solvol = solvol; %microliters
FAHvol = FAHvol; %microliters

%reference indices for the organic, and solvent. The antisolvent,
%inorganic,and FAH identities remain constant.
antiind = 8;
solvind = find(strcmp(solventord(1,:),solventid) == 1);
orgind = find(strcmp(organicord(1,:),orgid) == 1);
FAHind = 6;
inorgind = 7;


antimols = 0;
orgmols = orgconc * (solvol * 10^-6); %also small subtarranean mammal
inorgmols = inorgconc * (solvol * 10^-6); %moles
solventvolume = ((solvol * 10^-6) - ((orgmols * organicconv(orgind)) + (inorgmols * solventconv(inorgind)))); %liters
solvmols = solventvolume * (solventconv(solvind)^-1); % moles
FAHmols = (solventconv(FAHind)^-1) * (FAHvol * 10^-6); %moles


orgvol = (orgmols * organicconv(orgind)); %liters
inorgvol = (inorgmols * solventconv(inorgind)); %liters
end


