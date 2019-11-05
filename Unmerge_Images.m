% -------------------------------------------------------------------------
%  Name: Unmerge_Images.m
%  Version: 1.0
%  Environment: Matlab 2019a
%  Date: 23/08/2019
%  Authors: Mads Bergholt, Conor Horgan
% -------------------------------------------------------------------------
function [varargout] = Unmerge_Images(LSQ_Reg_Model, Component_Indices, varargin)

% Extract the least-squares regression model scores for each component
Scores = LSQ_Reg_Model(:,Component_Indices);

% Calculate indices for extracting images
tmp = 0;
indice(1) = 0;
for i = 2:nargin-1
   indice(i) = tmp + size(varargin{i-1},1)*size(varargin{i-1},2);
   tmp = indice(i);
end

% Extract data into each image
for i = 1:nargin-2
    data = reshape(Scores(indice(i)+1:indice(i+1),1), size(varargin{i},1), size(varargin{i},2));
    data = data/sum(sum(data));
    evalc(['xs' num2str(i) '= data']);
    varargout{i} = eval(['xs' num2str(i)]);
    data = [];
end

end
