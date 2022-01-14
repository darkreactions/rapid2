
function [concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,bincen,upperbound] = binspawnv2(spawnT,concdcm,concfah,concinorg,concorg,concsolv,molsdcm,molssolv,molsinorg,molsorg,molsFAH,h,h0,molfractdcm,molfractsolv,molfractFAH,molfractinorg,molfractorg,bins,dcdt,dsdt,ctotal,htotal,upperboundtotal,bincentotal,mtime,bincen,upperbound)
spawnT = spawnT;
boi = find(h > spawnT);

if isempty(boi) == 0
    for i = boi(1)
        %Coordinate Updates, centers and upper bound
        ncentL = bincentotal(i,:) - (htotal(i,:)./4);
        ncentH = bincentotal(i,:) + (htotal(i,:)./4);
        bincentotal(i,:) = ncentH;
        bincentotal = cat(1,bincentotal(1:i-1,:),ncentL,bincentotal(i:end,:));
        bincen = bincentotal(:,length(bincentotal(1,:)));

        nboundL = upperboundtotal(i,:) - (htotal(i,:)./2);
        nboundH = upperboundtotal(i,:);
        upperboundtotal(i,:) = nboundH;
        upperboundtotal = cat(1,upperboundtotal(1:i-1,:),nboundL,upperboundtotal(i:end,:));
        upperbound = upperboundtotal(:,length(upperboundtotal(1,:)));
        % Heights
        h(i) = h(i)/2;
        h = cat(1,h(1:i),h(i:end));
        h0 = cat(1,h0(1:i),h0(i:end));
        % Bins
        bins = length(h);
        % Mols of Species
        molsdcm(i) = molsdcm(i)/2;
        molsdcm = cat(1,molsdcm(1:i),molsdcm(i:end));
        
        molssolv(i) = molssolv(i)/2;
        molssolv = cat(1,molssolv(1:i),molssolv(i:end));
        
        molsinorg(i) = molsinorg(i)/2;
        molsinorg = cat(1,molsinorg(1:i),molsinorg(i:end));
        
        molsorg(i) = molsorg(i)/2;
        molsorg = cat(1,molsorg(1:i),molsorg(i:end));
        
        molsFAH(i) = molsFAH(i)/2;
        molsFAH = cat(1,molsFAH(1:i),molsFAH(i:end));
        
        % Concentrations of Species
        concdcm = cat(1,concdcm(1:i),concdcm(i:end));
        concfah = cat(1,concfah(1:i),concfah(i:end));
        concinorg = cat(1,concinorg(1:i),concinorg(i:end));
        concorg = cat(1,concorg(1:i),concorg(i:end));
        concsolv = cat(1,concsolv(1:i),concsolv(i:end));
        
        % Concentration Change Update
        molfractdcm = cat(1,molfractdcm(1:i),molfractdcm(i:end));
        molfractsolv = cat(1,molfractsolv(1:i),molfractsolv(i:end));
        molfractFAH = cat(1,molfractFAH(1:i),molfractFAH(i:end));
        molfractinorg = cat(1,molfractinorg(1:i),molfractinorg(i:end));
        molfractorg = cat(1,molfractorg(1:i),molfractorg(i:end));
        
        
        dcdt = zeros(bins,1);
        dsdt = zeros(bins,1);
        
        % Running Total Update
        ctotal = cat(1,ctotal(1:i,:),ctotal(i:end,:));
        
        htotal(i,:) = htotal(i,:)/2;
        htotal = cat(1,htotal(1:i,:),htotal(i:end,:));
        % Troubleshooting Display
        %disp(mtime)
        %disp(boi)
        %disp('Spawn!')
    end
end

%parameters = 'C:\Users\Rod Keesey\Documents\MATLAB\DiffusionModel\Tmp\binmergeparameters.mat';
%save(parameters,'spawnT','concdcm','concfah','concinorg','concorg','concsolv','molsdcm','molssolv','molsinorg','molsorg','molsFAH','h','h0','molfractdcm','molfractsolv','bins','dcdt','dsdt','ctotal','htotal','upperboundtotal','bincentotal','mtime','bincen','upperbound')

end

