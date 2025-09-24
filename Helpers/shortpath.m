function [x,y] = shortpath(x,y,DEM,iter)

    if nargin == 3, iter = 1.42; end

    xi = x; yi = y;

    % Concatenate x and y to form coordinate matrix 'c'
    c = [x(:), y(:)];
    
    % Find unique coordinates 'uc' without sorting
    uc = unique(c, 'rows', 'stable');
    
    % Extract unique x and y coordinates
    x = uc(:, 1);
    y = uc(:, 2);
    
    s = [];
    t = [];
    d = [];
    
    if size(x,1) == 1
        x = x';
        y = y';
    end
    
    runLoop = true;
    
    while runLoop && iter < 10
        try
            [P,D] = rangesearch([x y],[x y],DEM.cellsize*iter);
            
            for i3 = 1:numel(x)
                if numel(D{i3})>0
                    ud = unique(D{i3},'stable');
                    indi = P{i3}(D{i3}==ud(2));
                    di = D{i3}(D{i3}==ud(2));
            
                    if numel(indi) == 1 && numel(ud)>2, ind = [indi P{i3}(find(D{i3}==ud(3),1))]; 
                        d = cat(2,d,[di D{i3}(find(D{i3}==ud(3),1))]); 
                    else, ind = indi; d = cat(2,d,di); 
                    end
                   
                    s = cat(2,s,repmat(i3,1,numel(ind)));
                    t = cat(2,t,ind);
                end
            end

            G = graph(s,t,d);
    
            di = nan(1,numel(x));
            ind = nan(1,numel(x));
            for i3 = 1:numel(x)
                [~,di2] = shortestpathtree(G,i3); 
                if any(isinf(di2)), error(); end
                [maxi,indi] = max(di2);
                ind(i3) = indi;
                di(i3) = maxi;
            end

            [~,si] = max(di);
            ti = ind(si);
            
            %     P = flip(shortestpath(G,si,ti));
            P = shortestpath(G,si,ti);

            x = x(P);
            y = y(P);
    
            runLoop = false;
        catch
            iter = iter+.5;
        end
    end  
    
    if numel(x) == 0
        warning('Error during the shortpathing. Return the initial coordinates.')
        x = xi; y = yi;
    end

end
