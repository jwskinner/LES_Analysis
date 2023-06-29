function varargout = loadNetCDF(filename, varargin)
% LOADNETCDF Load specified variables from a netCDF file
%   [var1, var2, ..., varn] = loadNetCDF(filename, var1, var2, ..., varn) loads the
%   specified variables from the netCDF file `filename` and returns them as output arguments.
%
%   Example:
%   [lat, lon, temp] = loadNetCDF('data.nc', 'latitude', 'longitude', 'temperature');

finfo = ncinfo(filename); 

% Loop through each specified variable
for i = 1:nargin-1
    clear varargin{i} % Clear variable from memory if its held there already
    % Load the variable using ncread
    varargout{i}.Data = ncread(filename, varargin{i});
    varargout{i}.Units = ncreadatt(filename, varargin{i}, 'units');
    varargout{i}.Desc = ncreadatt(filename, varargin{i}, 'description');
    varargout{i}.Name = varargin{i};

end
