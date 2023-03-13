
function [XXdum, XYdum,YYdum ] = dsge_prior(DSGE_type, para,p);

% nshock = 9;

if DSGE_type ==1
   [ZZ] = make_zz_NK();
   [XXdum, XYdum, YYdum] = Make_dsge_moments_NK(para,ZZ,p);

elseif DSGE_type ==2
   [ZZ] = make_zz_FF();
   [XXdum, XYdum, YYdum] = Make_dsge_moments_FF(para,ZZ,p);
end


