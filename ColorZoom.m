function ColorZoom(zoom, varargin)
%ColorZoom(zoom)  Zoom in/out a color axis. >1 = produce larger color range
clim = get(gca,'clim');
crange = diff(clim);
clim(2) = clim(1) + crange .* zoom;
set(gca,'clim',clim);