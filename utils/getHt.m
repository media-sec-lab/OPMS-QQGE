function [ Ht] = getHt(P)  
        P = P(:);
        H = -((P).*log2(P));
        H((P<eps) | (P > 1-eps)) = 0;
        Ht = sum(H);
end

