function T = statsT(T)

%     n1 = numel(T.int.xi);
%     xii = cell(n1,1);
%     yii = cell(n1,1);
%     zii = cell(n1,1);
%     dii = cell(n1,1);
%     
%     for i1 = 1:n1
%         n2 = numel(T.int.xi{i1});
%         xii{i1} = cell(n2,1);
%         yii{i1} = cell(n2,1);
%         zii{i1} = cell(n2,1);
%         dii{i1} = cell(n2,1);
%     
%         for i2 = 1:n2
% %             mLen = max(cellfun('length', T.int.xi{i1}{i2}));
%             mLen = 0;
%             for j = 1:numel(T.int.xi{i1}{i2})
%                 mLen = min(mLen, length(T.int.xi{i1}{i2}{j}));
%             end
% 
%             xii{i1}{i2} = nan(1, mLen);
%             yii{i1}{i2} = nan(1, mLen);
%             zii{i1}{i2} = nan(1, mLen);
%             dii{i1}{i2} = nan(1, mLen);
%     
%             for idx = 1:mLen
% %                 vIdx = find(cellfun('length', T.int.xi{i1}{i2}) >= idx);
%                 vIdx = [];
%                 for j = 1:numel(T.int.xi{i1}{i2})
%                     if length(T.int.xi{i1}{i2}{j}) >= idx
%                         vIdx = [vIdx, j];
%                     end
%                 end
% 
% 
%                 if ~isempty(vIdx)
% %                     xiVals = cellfun(@(x) x(idx), T.int.xi{i1}{i2}(vIdx));
% %                     yiVals = cellfun(@(x) x(idx), T.int.yi{i1}{i2}(vIdx));
% %                     ziVals = cellfun(@(x) x(idx), T.int.zi{i1}{i2}(vIdx));
% %                     diVals = cellfun(@(x) x(idx), T.int.di{i1}{i2}(vIdx));
% %                     xiVals = zeros(1, numel(vIdx));
% %                     yiVals = zeros(1, numel(vIdx));
%                     ziVals = zeros(1, numel(vIdx));
%                     diVals = zeros(1, numel(vIdx));
%                     
%                     for k = 1:numel(vIdx)
%                         xiVals(k) = T.int.xi{i1}{i2}{vIdx(k)}(idx);
%                         yiVals(k) = T.int.yi{i1}{i2}{vIdx(k)}(idx);
%                         ziVals(k) = T.int.zi{i1}{i2}{vIdx(k)}(idx);
%                         diVals(k) = T.int.di{i1}{i2}{vIdx(k)}(idx);
%                     end
%    
%                     xii{i1}{i2}(idx) = xiVals(1);
%                     yii{i1}{i2}(idx) = yiVals(1);
%                     zii{i1}{i2}(idx) = mean(ziVals);
%                     dii{i1}{i2}(idx) = mean(diVals);
%                 end
%             end
%         end
%     end

%%
% for i1 = 1:n1
%     n2 = numel(T.int.xi{i1});
%     xii{i1} = cell(n2, 1);
%     yii{i1} = cell(n2, 1);
%     zii{i1} = cell(n2, 1);
%     dii{i1} = cell(n2, 1);
% 
%     for i2 = 1:n2
%         midx = ceil(numel(T.int.xi{i1}{i2}) / 2);
%         xi = T.int.xi{i1}{i2}{midx};
%         yi = T.int.yi{i1}{i2}{midx};
%         zi = T.int.zi{i1}{i2}{midx};
%         
%         xii{i1}{i2} = xi;
%         yii{i1}{i2} = yi;
%         dii{i1}{i2} = T.int.di{i1}{i2}{midx};
% 
%         % Collect all coordinates and values from the current cell
%         xi_all = vertcat(T.int.xi{i1}{i2}{:});
%         yi_all = vertcat(T.int.yi{i1}{i2}{:});
%         zi_all = vertcat(T.int.zi{i1}{i2}{:});
% 
%         coords_all = [xi_all, yi_all];
% 
%         % Initialize cell arrays to store the indices of the closest points
%         closest_indices = cell(numel(xi), 1);
% 
%         % Calculate distances from each point in coords_all to all xi, yi points
%         dists = pdist2(coords_all, [xi,yi]);
% 
%         % Find the closest xi, yi for each point in coords_all
%         [~, minIdx] = min(dists, [], 2);
% 
%         % Assign each element in coords_all to the closest xi, yi point
%         for j = 1:numel(minIdx)
%             closest_indices{minIdx(j)} = [closest_indices{minIdx(j)}; j];
%         end
% 
%         % Calculate mean values for zii and dii using the assigned points
%         zii{i1}{i2} = nan(size(xi));
% 
%         for k = 1:numel(closest_indices)
%             zii{i1}{i2}(k) = mean([zi(k);zi_all(closest_indices{k})],'omitnan');
%         end
%     end
% end

    %%

    xii = cell(2,1);
    yii = cell(2,1);
    zii = cell(2,1);
    dii = cell(2,1);

    for i1 = 1:2
        n2 = numel(T.int.x{i1});
        xii{i1} = cell(n2, 1);
        yii{i1} = cell(n2, 1);
        zii{i1} = cell(n2, 1);
        dii{i1} = cell(n2, 1);
    
        for i2 = 1:n2
            if ~isempty(T.int.x{i1}{i2})
                midx = ceil(numel(T.int.x{i1}{i2})/2);
                xii{i1}{i2} = T.int.x{i1}{i2}{midx};
                yii{i1}{i2} = T.int.y{i1}{i2}{midx};
                zii{i1}{i2} = T.int.z{i1}{i2}{midx};
                dii{i1}{i2} = T.int.d{i1}{i2}{midx};
            end
        end
    end

    for i2 = 1:numel(T.bl.x)

        try
            T.stats.x{i2,1} = [flip(xii{1}{i2});xii{2}{i2}(2:end)]; T.stats.y{i2,1} = [flip(yii{1}{i2});yii{2}{i2}(2:end)];
            T.stats.z{i2,1} = [flip(zii{1}{i2});zii{2}{i2}(2:end)]; T.stats.d{i2,1} = [flip(dii{1}{i2}*(-1));dii{2}{i2}(2:end)];
        catch
            % Set all T.stats fields to empty if any input data is empty
            T.stats.x{i2,1} = []; T.stats.y{i2,1} = [];
            T.stats.z{i2,1} = []; T.stats.d{i2,1} = [];
        end
    end
    
    for i1 = 1:numel(T.bl.x)

        try        
            id1 = find(T.stats.d{i1}<0); id2 = find(T.stats.d{i1}>0);
            Z1 = T.stats.z{i1}(id1); Z2 = T.stats.z{i1}(id2);
            D1 = T.stats.d{i1}(id1); D2 = T.stats.d{i1}(id2);
        catch
        end

        try
            T.stats.Z1(i1) = max(Z1);
            G1 = atan(diff(Z1) ./ diff(abs(D1(:)))) * (180/pi); 
            T.stats.G1(i1) = mean(G1);
        catch
            T.stats.Z1(i1) = nan;
            T.stats.G1(i1) = nan;
        end
             
        try
            T.stats.Z2(i1) = max(Z2);
            G2 = atan(diff(Z2) ./ diff(abs(D2(:)))) * (180/pi); 
            T.stats.G2(i1) = mean(G2);
        catch
            T.stats.Z2(i1) = nan;
            T.stats.G2(i1) = nan;
        end

        try
            T.stats.W(i1) = T.stats.d{i1}(1)*(-1) + T.stats.d{i1}(end);
            T.stats.H(i1) = mean([T.stats.Z1(i1) T.stats.Z2(i1)]) - min(T.stats.z{i1});
        catch
            T.stats.W(i1) = nan;
            T.stats.H(i1) = nan;
        end
    end

%     figure, imageschs(T.DEM,[],'colormap',[.7 .7 .7],'colorbar',false), hold on
%     for i1 = 1:numel(T.bl.x) 
%         plot(T.stats.x{i1},T.stats.y{i1},'k'), hold on, drawnow
%     end
%     plot(T.bl.x, T.bl.y, '.-', 'color', 'r', 'LineWidth', 4); hold on
% 
%     d = sqrt(diff(T.bl.x).^2 + diff(T.bl.y).^2);
%     d = [0; cumsum(d)];
%     
%     
%     figure;
%     hold on;
%     
%     % Plot the first part with transparency
%     for i1 = 1:numel(T.bl.x)
%         x = T.stats.d{i1};
%         y = movmean(T.stats.z{i1}, 3);
%         plot3(x, ones(1, numel(x)) * d(i1), y, 'Color', [0 0 0 0.7]); hold on, 
%     end


end