% -------------------------------------------------------------------------
%  Name: Merge_Images.m
%  Version: 1.0
%  Environment: Matlab 2019a
%  Date: 22/08/2019
%  Authors: Mads Bergholt, Conor Horgan
% -------------------------------------------------------------------------

function [Merged_Output] = Merge_Images(varargin)

% Reshape input data to be (Num_Spectra x Spectrum_Length)
for i = 1:nargin
    Reshaped_Data = reshape(varargin{i}, size(varargin{i}, 1) * size(varargin{i}, 2), size(varargin{i}, 3));
    varargout{i} = Reshaped_Data;
    Reshaped_Data = [];
end

% Combine reshaped data into a single matrix for analysis
Merged_Output = cell2mat([varargout']);
Merged_Output = double(Merged_Output);

end
