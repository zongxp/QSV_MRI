function compare_cd(data,cd_fit,roi,range,range_angle)

if ~exist('roi','var')
    roi=ones(size(cd_fit));
end

if ~exist('range','var')
    
  range=min_max([imag(data(:));real(data(:))]);
end
if ~exist('range_angle','var')
    range_angle=[-180,180];
end

range_angle=range_angle/180*pi;

data(roi==0)=0;
cd_fit(roi==0)=0;


 an_data=(angle(data)-range_angle(1))/(range_angle(2)-range_angle(1))*(range(2)-range(1))+range(1);
 an_fit=(angle(cd_fit)-range_angle(1))/(range_angle(2)-range_angle(1))*(range(2)-range(1))+range(1);
 

  %im=cat(3,real(data),real(cd_fit),imag(data),imag(cd_fit),an_data,an_fit);
  
  im=cat(3,real(data),real(cd_fit),imag(data),imag(cd_fit));
  
  
  im=repmat2(im,20);
  if ~exist('range','var')
      range=min_max(im(:));
  end
  subplot(1,2,1);imshow4(im,range,[2,2],2);
  
  subplot(1,2,2);
  clr=lines(2);
  hold off;
  myPlot(real(data(:)),real(cd_fit(:)),'o',clr(1,:));
  hold on;
  myPlot(imag(data(:)),imag(cd_fit(:)),'>',clr(2,:));
  
  xlabel('Measured');
  ylabel('Fit');
  mm=min_max([imag(data(:));real(data(:))]);
  hold on;plot(mm,mm,'k-');
  xlim(mm*1.2);
  ylim(mm*1.2);
  
  legend('Real','Imaginary');
  
  