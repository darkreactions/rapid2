
function sse = DiffusionSSE3(k)

    filename = 'C:\Users\Rod\OneDrive\Documents\Matlab Stuff\DiffusionModel\Tmp\parameters.mat';
    load(filename)
    
    k=k;
    kin = k(1);
    kout = k(2);

    [t,hprofile,cprofile,columnheightout] = DiffusionAllSpeciesV2(columnheight,radius,m0,bins,dt,tduration,d,mva,cs0,co0,kin,kout,antimols,solvmols,inorgmols,orgmols,FAHmols,htov,solventord,organicord,solventconv,organicconv,antiind,solvind,orgind,FAHind,inorgind,vtotal);    
    
    % Making an array of the observed timepoints in the model data
    modeldata = cat(2,transpose(columnheightout(1:end-1)),transpose(t));
    d2comp = zeros(length(solventd),1);
    
    
    for i = 1:length(solventd)
        d2comp(i) = modeldata(find(modeldata(:,2)==solventd(i,2)),1);
    end
    d2comp = cat(2,d2comp,solventd);

    % Calculating sum of squared differences
    sse = sum((d2comp(:,1) - d2comp(:,2)).^2);
    
    figure(1)

    plot(d2comp(:,3),d2comp(:,1))
    hold on
    plot(d2comp(:,3),d2comp(:,2))
    legend({'Model','Experimental'},'Location','southeast')
    title(['Solution Heights vs. Time for ' rxnid])
    xlabel('Seconds')
    ylabel('Height (cm)')
    hold off
    
end