function [t,htotal,ctotal,columnheight,totalmolsdcm,totalvolfrac,bincentotal,concdcm,concfah,concorg,concinorg,concsolv,h,volfracoi] = DiffusionAllSpeciesV2(columnheight,radius,molsdcm0,bins,dt,duration,d,mva,corg,cinorg,kin,kout,startingmolsdcm,totalmolssolvent,totalmolsinorg,totalmolsorg,totalmolsFAH,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind,vtotal)


    
    radius = radius; %radius of column in cm
    molsdcm0 = molsdcm0; %inital mols of antisolvent in bins
    bins = bins; %number of bins
    dt = dt; %timestep in seconds
    duration = duration; %length of experiment in seconds
    d = d; %diffusion constant
    dbase = d; %base diffusion constant for other species, we use DCM as it provides an upper bound
    mva = .0639; %molar volume of antisolvent (DCM) liters/mol
    kin = kin; %condensation rate
    kout = kout; %evaporation rate
    htov = htov; % = .002966; % cm/microliter
    columnheight = vtotal * 1000 * htov; %height of column in cm
    cooldown = 1;

    % Species identities and conversion arrays
    antiind = antiind;
    solvind = solvind;
    orgind = orgind;
    inorgind = inorgind;
    FAHind = FAHind;
    solventord = solventord; %holds the reference ids for solventconv
    solventconv = solventconv; %converts moles to liters
    organicord = organicord; %holds the reference ids for orgconv
    organicconv = organicconv; %converts moles to liters


    % Creating initial arrays 
    
    %NOTE THAT THE TOP BIN IS THE LAST IN THE ARRAY
    %BECAUSE OF INDEXING

    h = ones(bins,1) .* (columnheight/bins); %creates an array of bins with the starting height, h
    h0 = h;
    molsdcm = ones(bins,1) .* molsdcm0; %creates starting mols of antisolvent, m
    %molsdcm(20) = .005; %test parameter, leave commented out
    molssolv = ones(bins,1) .* (totalmolssolvent / bins);
    %molssolv(19) = .0005; %test parameter, leave commented out
    molsinorg = ones(bins,1) .* (totalmolsinorg / bins);
    molsorg = ones(bins,1) .* (totalmolsorg / bins);
    molsFAH = ones(bins,1) .* (totalmolsFAH / bins);
    
    %inital rates of concentration change for all species
    dcdt = zeros(bins,1); 
    dsdt = zeros(bins,1);
    didt = zeros(bins,1);
    dodt = zeros(bins,1);
    dfdt = zeros(bins,1);
    t = 0:dt:duration; %creates array of timesteps

    %Height and Concentration Update
    h = h0 + (molsdcm.* (mva * 10^6 * htov)); %cm
    
    %Concentration from intital conditions
    concdcm = molsdcm ./ (h * htov^-1 / 10^6); %Updating Concentration of dcm
    concfah = molsFAH ./ (h * htov^-1 / 10^6);
    concinorg = molsinorg ./ (h * htov^-1 / 10^6);
    concorg = molsorg ./ (h * htov^-1 / 10^6);
    concsolv = molssolv ./ (h * htov^-1 / 10^6);
    
    %Updating Mol Fraction of Species
    molfractdcm = molsdcm ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH); %molfraction
    molfractsolv = molssolv ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    molfractFAH = molsFAH ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    molfractinorg = molsinorg ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    molfractorg = molsorg ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
    
    %Height and Concentration Running total Array
    htotal = []; %cm
    ctotal = [];% mols/liter
    
    %bin center initial
    bincentotal = []; %cm coordinate
    
    %upper bound initial
    upperboundtotal = []; %cm
    

    % Conservation test stuff
    columnheight = sum(h);
    totalmolsdcm = sum(molsdcm);
    molarity = sum(concdcm);
    bincheck = bins;
    totalvolfrac = [];
    intertotalvolfrac = [];
    %dmdttotal = [];
    poi = .3;
    toi = 1;
    spawnT = h0(1)*2;
 %%
    for mtime = 1:length(t)
        
        %for DCM 
        [molsdcm,dcdt,dmdt]= diffcalcV1(bins,d,h,concdcm,htov,kin,kout,molfractdcm,molsdcm,dt);
        %for Solvent
        [molssolv,dsdt,dmsdt]= diffcalcV1(bins,dbase,h,concsolv,htov,0,0,molfractsolv,molssolv,dt);
        %for FAH
        [molsFAH,dfdt,dmfdt]= diffcalcV1(bins,dbase,h,concfah,htov,0,0,molfractFAH,molsFAH,dt);
        %for inorganic solute
        [molsinorg,didt,dmidt]= diffcalcV1(bins,dbase,h,concinorg,htov,0,0,molfractinorg,molsinorg,dt);
        %for organic solute
        [molsorg,dodt,dmodt]= diffcalcV1(bins,dbase,h,concorg,htov,0,0,molfractorg,molsorg,dt);
        
        
        %Updating Heights
        hspecies = (molsFAH * solventconv(FAHind) + molsinorg * solventconv(inorgind)+ molsorg * organicconv(orgind) + molssolv * solventconv(solvind)).*(10^6 * htov);
        h = hspecies + (molsdcm.* (mva * 10^6 * htov)); %Calculating new height from mols
        
        
        
        %Species Concentration updates in Mols/Liter
        concdcm = molsdcm ./ (h * htov^-1 / 10^6); %Updating Concentration of dcm
        concfah = molsFAH ./ (h * htov^-1 / 10^6);
        concinorg = molsinorg ./ (h * htov^-1 / 10^6);
        concorg = molsorg ./ (h * htov^-1 / 10^6);
        concsolv = molssolv ./ (h * htov^-1 / 10^6);
        
        %Updating Molar Fraction of species
        molfractdcm = molsdcm ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
        molfractsolv = molssolv ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
        molfractFAH = molsFAH ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
        molfractinorg = molsinorg ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
        molfractorg = molsorg ./ (molsdcm + molssolv + molsinorg + molsorg + molsFAH);
        
        %Calculating Bin Centers
        bincen = h ./ 2; %% Finds 1 dimensional coordinate of each bin CENTER
        for t2height = 2:length(bincen)
            bincen(t2height) = (sum(h(1:t2height-1)) + bincen(t2height));
        end
        
        bincentotal = cat(2,bincentotal,bincen);
        
        %Calculating Bin Upper Boundary
        upperbound = h;
        for uppercalc = 2:length(upperbound)
            upperbound(uppercalc) = sum(h(1:uppercalc));
        end
        
        upperboundtotal = cat(2,upperboundtotal,upperbound);

        % Conservation Calculations
        ctotal = cat(2,ctotal,concdcm);
        htotal = cat(2,htotal,h);
        columnheight = cat(2,columnheight,sum(h)); %gives total height
        totalmolsdcm = cat(2,totalmolsdcm,sum(molsdcm)); %gives total mols
        molarity = cat(2,molarity,sum(concdcm));
        %dmdttotal = cat(2,dmdttotal,dmdt); %no bin spawning
        
        % Volume Fraction Calculations in bin,
        dcminbin = concdcm * solventconv(antiind);
        fahinbin = concfah * solventconv(FAHind);
        inorginbin = concinorg * solventconv(inorgind);
        orginbin = concorg * organicconv(orgind);
        solvinbin = concsolv * solventconv(solvind);

        volfrac = cat(2,dcminbin,fahinbin,inorginbin,orginbin,solvinbin, dcminbin + fahinbin + inorginbin + orginbin + solvinbin );
        totalvolfrac = cat(2,totalvolfrac,volfrac(length(bincen),:).');
        
        % Linear Interpolation and Volume Fraction
        [dcminterpolated,fahinterpolated,inorginterpolated,orginterpolated,solvinterpolated] = LinearInterp(bincen,concdcm,concfah,concinorg,concorg,concsolv,poi);
        dcminterfrac = dcminterpolated * solventconv(antiind);
        fahinterfrac = fahinterpolated * solventconv(FAHind);
        inorginterfrac = inorginterpolated * solventconv(inorgind);
        orginterfrac = orginterpolated * organicconv(orgind);
        solvinterfrac = solvinterpolated * solventconv(solvind);
        
        intervolfrac = cat(2,dcminterfrac,fahinterfrac,inorginterfrac,orginterfrac,solvinterfrac,dcminterfrac+fahinterfrac+inorginterfrac+orginterfrac+solvinterfrac);
        intertotalvolfrac = cat(2,intertotalvolfrac,intervolfrac.');
        
        %Bin Spawning
        if length(concdcm) < 20 && cooldown < 0 %elapsed < .005
            [concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,bincen,upperbound] = binspawnv2(spawnT,concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,mtime,bincen,upperbound);
            [concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,bincen,upperbound] = binspawnv2(spawnT,concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,mtime,bincen,upperbound);
            cooldown = 0; %1;
            %disp("merged")
            %disp(elapsed)
        end
        
        %Bin Merging
        if bins > 1
            [concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,bincen,upperbound] = binmergeV1(spawnT,concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,mtime,bincen,upperbound);  
        end
        
        cooldown = cooldown -1;
        %disp(cooldown)
        % Graphing
        %if mtime == mtime
        %if mtime == length(t)
        %    figure(1)
        %    subplot(3,2,1)
        %    scatter(bincen,concdcm)
        %    for line = 1:length(upperbound)
        %        xline(upperbound(line))
        %    end
        %    axis([0 sum(h) 0 inf])
        %    title("Concentration (mols/liter) per Height(cm)")
        %    subplot(3,2,2)
        %    plot(totalmolsdcm)
        %    axis([1 inf 0 inf])
        %    title("Total Moles of DCM per Timestep")
        %    subplot(3,2,3)
        %    plot(columnheight)
        %    axis([1 inf 0 inf])
        %    title("Total Height (cm) per Timestep")
        %    subplot(3,2,4)
        %    plot(totalvolfrac(1,:))
        %    hold on
        %    plot(totalvolfrac(2,:))
        %    plot(totalvolfrac(3,:))
        %    plot(totalvolfrac(4,:))
        %    plot(totalvolfrac(5,:))
        %    plot(totalvolfrac(6,:))
        %    axis([1 inf 0 1.5])
        %    hold off
        %    title("Bin Volume-Fraction of all Species per timestep")
        %    subplot(3,2,5)
        %    plot(intertotalvolfrac(1,:))
        %    hold on
        %    plot(intertotalvolfrac(2,:))
        %    plot(intertotalvolfrac(3,:))
        %    plot(intertotalvolfrac(4,:))
        %    plot(intertotalvolfrac(5,:))
        %    plot(intertotalvolfrac(6,:))
        %    axis([1 inf -2 2])
        %    hold off
        %    title("Interpolated Volume-Fraction of all Species per timestep")
            %subplot(3,2,6)
            %plot(1)
            %axis([0 4 0 inf])
            %rectangle('Position', [.99 0 2.02 2],'LineWidth',1)
            %rectangle('Position', [1 0 2 h(1,1)],'FaceColor',[0+dcminbin(1) 1-dcminbin(1) 1-dcminbin(1)],'EdgeColor',[0+dcminbin(1) 1-dcminbin(1) 1-dcminbin(1)])
            %for plotloop = 2:length(h(:,1))
            %    rectangle('Position',[1 sum(h(1:plotloop-1)) 2 h(plotloop)],'FaceColor',[0+dcminbin(plotloop) 1-dcminbin(plotloop) 1-dcminbin(plotloop)],'EdgeColor',[0+dcminbin(plotloop) 1-dcminbin(plotloop) 1-dcminbin(plotloop)])
            %    axis([0 4 0 inf])
            %end
            %title("Column Visualization of Model for DCM Concentration")
        %end
        %if mtime == mtime
        %if mtime == length(t)
        %    figure(2)
        %    subplot(3,2,1)
        %    plot(bincen,concdcm)
        %    for line = 1:length(upperbound)
        %        xline(upperbound(line))
        %    end
        %    axis([0 sum(h) 0 inf])
        %    title("DCM Concentration (mols/liter) per Height(cm)")
        %    subplot(3,2,2)
        %    plot(bincen,concsolv)
        %    for line = 1:length(upperbound)
        %        xline(upperbound(line))
        %    end
        %    axis([0 sum(h) 0 inf])
        %    title("Solvent Concentration (mols/liter) per Height(cm)")
        %    subplot(3,2,3)
        %    plot(bincen,concfah)
        %    for line = 1:length(upperbound)
        %        xline(upperbound(line))
        %    end
        %    axis([0 sum(h) 0 inf])
        %    title("FAH Concentration (mols/liter) per Height(cm)")
        %    subplot(3,2,4)
        %    plot(bincen,concinorg)
        %    for line = 1:length(upperbound)
        %        xline(upperbound(line))
        %    end
        %    axis([0 sum(h) 0 inf])
        %    title("Inorganic Concentration (mols/liter) per Height(cm)")
        %    subplot(3,2,5)
        %    plot(bincen,concorg)
        %    for line = 1:length(upperbound)
        %        xline(upperbound(line))
        %    end
        %    axis([0 sum(h) 0 inf])
        %    title("Organic Concentration (mols/liter) per Height(cm)")
        %end
        %if mtime == mtime
        %if mtime == length(t)
        %    figure(3)
        %    plot(1)
        %    axis([0 4 0 inf])
        %    rectangle('Position', [.99 0 2.02 2],'LineWidth',1)
        %    rectangle('Position', [1 0 2 h(1,1)],'FaceColor',[0+dcminbin(1) 1-dcminbin(1) 1-dcminbin(1)],'EdgeColor',[0+dcminbin(1) 1-dcminbin(1) 1-dcminbin(1)])
        %    for plotloop = 2:length(h(:,1))
        %        rectangle('Position',[1 sum(h(1:plotloop-1)) 2 h(plotloop)],'FaceColor',[0+dcminbin(plotloop) 1-dcminbin(plotloop) 1-dcminbin(plotloop)],'EdgeColor',[0+dcminbin(plotloop) 1-dcminbin(plotloop) 1-dcminbin(plotloop)])
        %        axis([0 4 0 inf])
        %    end
        %    title("Column Visualization of Model for DCM Concentration")
        %end
    end
    volfracoi = intertotalvolfrac(:,toi);
end