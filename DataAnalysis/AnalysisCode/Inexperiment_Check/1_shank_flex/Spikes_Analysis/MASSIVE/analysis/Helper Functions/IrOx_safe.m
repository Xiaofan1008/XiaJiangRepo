%% Calculates the safe charge injection limit for IrOx
function safe_charge = IrOx_safe(diameter)
% diameter is expected in units of um
area = pi*(diameter/2)^2;
area_cm = area*1e-8;
safe_charge = (area_cm * 10^1.5)^0.5;
end