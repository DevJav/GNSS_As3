function [d] = saastamoinensaastamoinen(elev, lambda, height, P, T, RH)
    e_s = 6.106 * RH * exp((17.15*T-4684)/(T-38.46));
    D = 1 + 0.0026 * cos(2*height) + 0.00028 * lambda;
    d_zenith = 0.002277 * D * (P + (1255/T+0.05) * e_s);
    d = 1 / sin(elev) * d_zenith;
end