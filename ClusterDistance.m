function r = ClusterDistance(C, varargin)
%convert Cluster.xy to radius/angle
%

if isfield(C,'cx') && isfield(C,'cy')
        rx = C.xyr(3);
        ry = C.xyr(4);
        if isfield(C,'xy')
            xy = C.xy;
        else
            xy(:,1) = C.cx;
            xy(:,2) = C.cy;
        end
        if isfield(C,'aspectratio') & C.aspectratio > 0
            xys = xyrotate(C.cx-C.xyr(1),(xy(:,2)-C.cy) ./C.aspectratio,C.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./C.aspectratio)).^2;
        else
            xys = xyrotate(xy(:,1)-C.xyr(1),xy(:,2)-C.xyr(2),C.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;
        end
        if nargout > 1
            y = xys(:,2);
        end
elseif isfield(C,'xy')
        rx = C.xyr(3);
        ry = C.xyr(4);
        xy = C.xy;
        if isfield(C,'aspectratio') & C.aspectratio > 0
            xys = xyrotate(xy(:,1)-C.xyr(1),(xy(:,2)-C.xyr(2)) ./C.aspectratio,C.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./C.aspectratio)).^2;
        else
            xys = xyrotate(xy(:,1)-C.xyr(1),xy(:,2)-C.xyr(2),C.angle);
            r = ((xys(:,1))./rx).^2 + ((xys(:,2))./ry).^2;
        end
        if nargout > 1
            y = xys(:,2);
        end
else
end