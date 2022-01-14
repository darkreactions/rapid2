function [solventmolarity,FAHmolarity,molsolv,molsFAH,vtotal] = StartingCalculationsV2(solventtypein,solvolin,FAHvolin,vi,cinorg,corg,amine)

%starting conditions
gblden = 1.12; %g/ml
dmfden = .948;
dmsoden = 1.1;
gbldmfden = 1.034;
dmfdmsoden = 1.024;
fahden = 1.22;
ki2den = 6.16;

%molar masses
gblmm = 86.090; %g/mol
dmfmm = 73.095;
dmsomm = 78.13;
gbldmfmm = 79.593;
dmfdmsomm = 75.613;
fahmm = 46.03;
ki2mm = 461.01;

%Amine Density
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

%Amine Molar Mass
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


%Values Passed through
solventtype = solventtypein;
solvol = solvolin / 1000; % given in microliters converted to ml
FAHvol = FAHvolin /1000; % given in microliters converted to ml
vi = vi /1000; %microliters covnerted to ml
cinorg = cinorg; %mols/L
corg = corg; %mols/L
amine = amine;




atarray = [gblden,dmfden,dmsoden,gbldmfden,dmfdmsoden,fahden,ki2den;gblmm,dmfmm,dmsomm,gbldmfmm,dmfdmsomm,fahmm,ki2mm];
ord = {'GBL','DMF','DMSO','GBL:DMF','DMF:DMSO','FAH','KI2'};

amarray = [Pden,Cyden,Acden,Etden,Meden,Phden,Prden,NNeden,NNpden,Beden;Pmm,Cymm,Acmm,Etmm,Memm,Phmm,Prmm,NNemm,NNpmm,Bemm];
amord = {'2Pyrrolidin1ium1ylethylammoniumiodide','CyclohexylmethylammoniumIodide','AcNH3I','EtNH3I','MeNH3I','PhEtNH3I','Propane13diammoniumIodide','NNDimethylethane12diammoniumiodide','NNDiethylpropane13diammoniumiodide','Benzenediaminedihydroiodide'};

%%
%Calculating Volumes


%Calculating Concentration of Solvent and FAH
mdarray = atarray(1,:) ./ atarray(2,:);
mdoi = mdarray(find(strcmp(ord(1,:),solventtype) == 1));
inmdam = amarray(2,:) ./ amarray(1,:); % ml organic /mol

%lx = ki2mm/ki2den; % ml inorganic/mol
%ly =  inmdam(find(strcmp(amord(1,:),amine) == 1));

%syms x y
%[solx,soly] = solve( x == (vi + x*lx + y * ly)*cinorg/1000, y == (vi + x*lx + y * ly)*corg/1000 );

%molsinorg = double(solx);
%molsorg = double(soly);

%volinorg = molsinorg * lx; %ml
%volorg = molsorg * ly; %ml

%vtotal = vi + volinorg + volorg;

vtotal = vi;

molsolv = solvol * mdoi; %Initial Mols of solvent
molsFAH = FAHvol * mdarray(6); %Initial Mols of FAH
%concentrations in Mols/Liter
solventmolarity = molsolv / vtotal *1000;
FAHmolarity = molsFAH / vtotal * 1000;

end
