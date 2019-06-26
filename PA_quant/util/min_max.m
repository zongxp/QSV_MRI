function a = min_max(data)
%a = min_max(data) 
% a is a 1 by 2 row vector containing min and max values of data.
 tmp = data(:);
 a = zeros(1,2);
 a(1) = min(tmp);
 a(2) = max(tmp);