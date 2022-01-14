

function [molsspecies,dcdtspecies,dmdtspecies]= diffcalcV1(bins,d,h,concspecies,htov,kin,kout,molfractspecies,molsspecies,dt)

%disp(bins)

if bins > 1
    for mpos = 2:bins-1
        dcdtspecies(mpos) = d/h(mpos) * ((concspecies(mpos+1)-concspecies(mpos))/(h(mpos+1)/2+(h(mpos)/2))) - (d/h(mpos) * ((concspecies(mpos)-concspecies(mpos-1))/(h(mpos)/2+(h(mpos-1)/2))));
    end
    dcdtspecies(bins) =(d/h(bins) * ((concspecies(bins)-concspecies(bins-1))/(h(bins)/2+(h(bins-1)/2)))); %top bin

    dcdtspecies(1) = d/h(1) * ((concspecies(2)-concspecies(1))/(h(2)/2+(h(1)/2))); %bottom bin
else
   dcdtspecies(:) = 0;
end


% Better Method of calculating change of mols per timestep, CONFIRMED WORKING
dmdtspecies = (h * htov^-1 / 10^6) .* dcdtspecies'; %all bins but top
dmdtspecies(bins) = kin - dmdtspecies(bins) -kout * molfractspecies(bins); %top bin

%Updating Running Totals of Mols
molsspecies = molsspecies + dmdtspecies * dt;



for checkneg = 1:length(molsspecies)
    if molsspecies(checkneg) < 0
        molsspecies(checkneg) = 0;
        %disp("Mass may not be conserved!")
    end
end

