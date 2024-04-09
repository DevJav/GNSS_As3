function [eph, ion] = brdc2eph(brdc,t0)

eph = struct;
ion = struct;

%BRDC2EPH Summary of this function goes here
%   Detailed explanation goes here
    for i=1:length(brdc.sat)

        [val,idx] = min(abs(t0 - brdc.eph{1,i}(23,:)));

        str = sprintf('PRN%02d',brdc.sat(i));
        eph.(str).t_oe     = brdc.eph{1,i}(23,idx);
        eph.(str).sqrtA    = sqrt(brdc.eph{1,i}(22,idx));
        eph.(str).e        = brdc.eph{1,i}(20,idx);
        eph.(str).deltan  = brdc.eph{1,i}(17,idx);
        eph.(str).M_0       = brdc.eph{1,i}(18,idx);
        eph.(str).omega_0    = brdc.eph{1,i}(25,idx);
        eph.(str).omegaDot   = brdc.eph{1,i}(30,idx);
        eph.(str).omega    = brdc.eph{1,i}(29,idx);
        eph.(str).i_0       = brdc.eph{1,i}(27,idx);
        eph.(str).iDot       = brdc.eph{1,i}(31,idx);
        eph.(str).C_rs      = brdc.eph{1,i}(16,idx); 
        eph.(str).C_rc      = brdc.eph{1,i}(28,idx); 
        eph.(str).C_is      = brdc.eph{1,i}(26,idx); 
        eph.(str).C_ic      = brdc.eph{1,i}(24,idx); 
        eph.(str).C_us      = brdc.eph{1,i}(21,idx); 
        eph.(str).C_uc      = brdc.eph{1,i}(19,idx); 
        eph.(str).a_f0      = brdc.eph{1,i}(12,idx);
        eph.(str).a_f1      = brdc.eph{1,i}(13,idx);
        eph.(str).a_f2      = brdc.eph{1,i}(14,idx);
        eph.(str).T_GD      = brdc.eph{1,i}(37,idx);
        eph.(str).t_oc      = brdc.eph{1,i}(39,idx);


    end

    ion.alpha0 = brdc.hdr.ionoAlpha(1);
    ion.alpha1 = brdc.hdr.ionoAlpha(2);
    ion.alpha2 = brdc.hdr.ionoAlpha(3);
    ion.alpha3 = brdc.hdr.ionoAlpha(4);

    ion.beta0 = brdc.hdr.ionoBeta(1);
    ion.beta1 = brdc.hdr.ionoBeta(2);
    ion.beta2 = brdc.hdr.ionoBeta(3);
    ion.beta3 = brdc.hdr.ionoBeta(4);

end

