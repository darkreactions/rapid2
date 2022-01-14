function [solvinterpolated] = LinearInterpV2(bincen,concsolv,poi)
bincen = bincen;
concsolv = concsolv;
poi = poi;


%% Finding Bins of Interest
if length(bincen) > 1
    boi1 = find(bincen <= poi, 1, 'last' );%% finding bin of interest, lower
    boi2 = find(bincen <= poi, 1, 'last' ) + 1; %% finding bin of interest, upper

    if poi < bincen(1)
        solvinterpolated = concsolv(1,:);
    elseif boi2 > length(bincen)
        solvinterpolated = concsolv(length(bincen),:);
    else
        solvinterpolated = diffinterp(boi1,boi2,bincen,concsolv,poi);
    end
else
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