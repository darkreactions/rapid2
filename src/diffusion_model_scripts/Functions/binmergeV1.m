
function [concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,bincen,upperbound] = binmergeV1(spawnT,concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,mtime,bincen,upperbound)

htov = .002966; % cm/microliter
spawnT = spawnT;
boi = find(h < h0/3);

if isempty(boi) == 0

    tomerge = boi(1);

    if tomerge == 1
        lower = 1;
        upper = 2;
    elseif tomerge == length(h)
        lower = length(h)-1;
        upper = length(h);
    else
        testvals = [tomerge+1,tomerge-1;h(tomerge+1),h(tomerge-1)];
        [lowvalue,index] = min(testvals(2,:));
        mergeto = testvals(1,index);
        upper = max([tomerge,mergeto]);
        lower = min([tomerge,mergeto]);
    end


    %% coordinate update (bincenter)

    bincentotalnew = upperboundtotal(upper,:) - ((htotal(lower,:) + htotal(upper,:))/2);
    bincentotal = cat(1,bincentotal(1:lower-1,:),bincentotalnew,bincentotal(upper+1:end,:));
    bincen = bincentotal(:,length(bincentotal(1,:)));

    %% coordinate update (upperbound)

    upperboundtotalnew = upperboundtotal(upper,:);
    upperboundtotal = cat(1,upperboundtotal(1:lower-1,:),upperboundtotalnew,upperboundtotal(upper+1:end,:));
    upperbound = upperboundtotal(:,length(upperboundtotal(1,:)));

    %% Height update
    htotalold = htotal;
    htotal = cat(1,htotal(1:lower-1,:),htotal(lower,:)+htotal(upper,:),htotal(upper+1:end,:));
    h = htotal(:,length(htotal(1,:)));
    h0 = cat(1,h0(1:lower-1,:),h0(lower,:)+h0(upper,:),h0(upper+1:end,:));
    %% bins update

    bins = length(h);

    %% Mols of species update

    molsdcm = cat(1,molsdcm(1:lower-1,:),molsdcm(lower,:)+molsdcm(upper,:),molsdcm(upper+1:end,:));
    molssolv = cat(1,molssolv(1:lower-1,:),molssolv(lower,:)+molssolv(upper,:),molssolv(upper+1:end,:));
    molsinorg = cat(1,molsinorg(1:lower-1,:),molsinorg(lower,:)+molsinorg(upper,:),molsinorg(upper+1:end,:));
    molsorg = cat(1,molsorg(1:lower-1,:),molsorg(lower,:)+molsorg(upper,:),molsorg(upper+1:end,:));
    molsFAH = cat(1,molsFAH(1:lower-1,:),molsFAH(lower,:)+molsFAH(upper,:),molsFAH(upper+1:end,:));

    %% Concentration of Species update

    concdcm = molsdcm ./ (h * htov^-1 / 10^6); %Updating Concentration of dcm
    concfah = molsFAH ./ (h * htov^-1 / 10^6);
    concinorg = molsinorg ./ (h * htov^-1 / 10^6);
    concorg = molsorg ./ (h * htov^-1 / 10^6);
    concsolv = molssolv ./ (h * htov^-1 / 10^6);

    volumetotal = (htotalold * htov^-1 / 10^6);

    newconc = ((volumetotal(lower,:) .* ctotal(lower,:)) + (volumetotal(upper,:) .* ctotal(upper,:)))./(htotal(lower,:) * htov^-1 / 10^6);
    ctotal = cat(1,ctotal(1:lower-1,:),newconc,ctotal(upper+1:end,:));

    %% Mole Fraction Changes
    molfractdcm = molsdcm ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH); %molfraction
    molfractsolv = molssolv ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    molfractFAH = molsFAH ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    molfractinorg = molsinorg ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    molfractorg = molsorg ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    %% DCDT Update
    dcdt = zeros(bins,1);
    dsdt = zeros(bins,1);

    %% Troubleshooting
    %disp(mtime)
    %disp(boi)
    %disp('Merge!')
end