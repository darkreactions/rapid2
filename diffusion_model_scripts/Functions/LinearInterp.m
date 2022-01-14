
function [dcminterpolated,fahinterpolated,inorginterpolated,orginterpolated,solvinterpolated] = LinearInterp(bincen,concdcm,concfah,concinorg,concorg,concsolv,poi)
bincen = bincen;
concdcm = concdcm;
concfah = concfah;
concinorg = concinorg;
concorg = concorg;
concsolv = concsolv;

poi = poi;


%% Finding Bins of Interest
if length(bincen) > 1
    boi1 = find(bincen <= poi, 1, 'last' );%% finding bin of interest, lower
    boi2 = find(bincen <= poi, 1, 'last' ) + 1; %% finding bin of interest, upper

    if poi < bincen(1)
        boi1 = 1;
        boi2 = 2;
    end

    if boi2 > length(bincen)
        boi1 = length(bincen)-1;
        boi2 = length(bincen);
    end

    dcminterpolated = diffinterp(boi1,boi2,bincen,concdcm,poi);
    fahinterpolated = diffinterp(boi1,boi2,bincen,concfah,poi);
    inorginterpolated = diffinterp(boi1,boi2,bincen,concinorg,poi);
    orginterpolated = diffinterp(boi1,boi2,bincen,concorg,poi);
    solvinterpolated = diffinterp(boi1,boi2,bincen,concsolv,poi);
else
    dcminterpolated = concdcm;
    fahinterpolated = concfah;
    inorginterpolated = concinorg;
    orginterpolated = concorg;
    solvinterpolated = concsolv;
end
%% Calculating the interpolated concentrations
function [concentrationout] = diffinterp(boi1,boi2,bincen,conc,poi)

m = (conc(boi2) - conc(boi1)) / (bincen(boi2) - bincen(boi1));
y_i = conc(boi2);
x_i = bincen(boi2);
b = y_i - (m * x_i);
concentrationout = m * poi + b;

end
end

