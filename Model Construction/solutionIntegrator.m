function [P] = solutionIntegrator(admsbl_p)
    
    [~,~,K] = size(admsbl_p);
    P = double.empty;
    for i = 1 : K
        new = admsbl_p(find(all(admsbl_p(:,:,i),2)),:,i);
        if ~isempty(P)
            I = find(~ismember(new,P,'rows'));
            P = [P;admsbl_p(I,:,i)];
        else
            P = new;
        end
    end

end

