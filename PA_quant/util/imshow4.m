function img=imshow4(img,range,shape,gap)
% imshow4(img, [ range, [shape )
%
% function to display series of images as a montage.
% 
% img - a 3D array representing a series of images
% range - window level or color map (similar to imshow)
% shape - a 2x1 vector representing the shape of the motage
%
% Example: 
% 		im = repmat(phantom(128),[1,1,6]);
%		figure;
%		imshow3(im,[],[2,3]);
%
% (c) Michael Lustig 2012

if ~exist('gap','var')
    gap=0;
end

img = img(:,:,:);
[sx,sy,nc] = size(img);

if nargin < 2 
	range = double([min(img(:)), max(img(:))]);
end

if isempty(range)
	range = [min(img(:)), max(img(:))];
    
    
end
if size(range,2)==3
  range=double(range);

  range=cat(1,range,[1,1,1]);
else
    img(img>range(2))=range(2);
end

  img=cat(1,img,ones(gap,sy,nc)*(1+range(2)));
  img=cat(2,img,ones(sx+gap,gap,nc)*(1+range(2)));

  
sx=sx+gap;
sy=sy+gap;

if nargin < 3

	if  ceil(sqrt(nc))^2 ~= nc;
   	   nc = ceil(sqrt(nc))^2;
   	   img(end,end,nc)=0;
	end


	img = reshape(img,sx,sy*nc);
	img = permute(img,[2,3,1]);
	img = reshape(img,sy*sqrt(nc),sqrt(nc),sx);
	img = permute(img,[3,2,1]);
	img = reshape(img,sx*sqrt(nc),sy*sqrt(nc));

else
    if nc<shape(1)*shape(2)
        img(:,:,end+1:shape(1)*shape(2))=0;
        nc=shape(1)*shape(2);
    end
	img = reshape(img,sx,sy*nc);
	img = permute(img,[2,3,1]);
	img = reshape(img,sy*shape(2),shape(1),sx);
	img = permute(img,[3,2,1]);
	img = reshape(img,sx*shape(1),sy*shape(2));
end

%imagesc(img,range); colormap(gray(256));axis('equal');
if nargout==0
imshow(img(1:end-gap,1:end-gap),range);
else
   res= img(1:end-gap,1:end-gap);
    
end
