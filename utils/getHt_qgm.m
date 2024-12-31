function [Ht] = getHt_qgm(P)   %in nats
    P = P(:);
    H = -((P).*log(P));
    H((P<eps) | (P > 1-eps)) = 0;
    Ht = sum(H);
end