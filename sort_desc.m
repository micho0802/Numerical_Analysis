function [res,index] = sort_desc(data);
% sort_desc Sort in descending order.
% res = sort_desc(data) sorts the elements of the vector
% data in descending order. [res,index] = sort(data) also
% returns an index vector index, i.e., res = data(index).
% When data is complex, the elements are sorted by
% abs(data).
% check usage
error(nargchk(1,1,nargin));
% sort data in ascending order using the MATLAB built-in
% function sort
[res,index] = sort(data);
% reorder the data
[m,n] = size(data);
if (m <= n)
  res = fliplr(res);
  index = fliplr(index);
else
  res = flipud(res);
  index = flipud(index);
end