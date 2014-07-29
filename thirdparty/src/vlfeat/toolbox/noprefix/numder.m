function varargout = numder(varargin)
% VL_NUMDER  Numerical derivative
%   D = VL_NUMDER(FUNC, X) computes the numerical derivative of 
%   the function FUNC at point X.
%
%   D = VL_NUMDER(FUNC, X, ARG1, ARG2, ...) allow to pass extra
%   parameters to the function FUNC.
%
%   See also:: VL_NUMDER2(), VL_HELP().
[varargout{1:nargout}] = vl_numder(varargin{:});
