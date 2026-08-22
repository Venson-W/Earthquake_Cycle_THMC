function zeta = compute_zeta_V(tauf,eta,V)

    s = eta .* V;

    zeta = tauf + s;
end