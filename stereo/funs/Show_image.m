function Show_image(img,sz)
%show_image(img,sz)
%
%Shows an image of size 'sz'.  By default, 'sz' is the true size
%of the image.

if nargin==1
	sz=size(img);
end

imagesc(img);
colormap(gray);
axis off
%truesize(sz);
