function [p] = searchSpace(CG,s,p,org_sample)

    D = pdist2(CG,p);
    [~,I] = sort(D);
    p_sorted = p(I',:);
    [r,~] = size(p);
    if r > 0.1 * org_sample
        exploit = p_sorted(1:round(0.1 * org_sample),:);
        II = I(1:round(0.1 * org_sample));
        p(II',:) = [];
        D = pdist2(s,exploit);
        [~,I] = sort(D);
        exploit = exploit(I',:);
        p = [exploit;p];
    else
        exploit = p_sorted;
        D = pdist2(s,exploit);
        [~,I] = sort(D);
        exploit = exploit(I',:);
        p = exploit;
    end
      
end

