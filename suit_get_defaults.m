function varargout = suit_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = suit_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT spm_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit spm_defaults.m.
% v.2.5: Compatibility with SPM job manager (26/08/2010) 
%_______________________________________________________________________
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)

global defaults_suit;
if isempty(defaults_suit)
    suit_defaults;
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');
if nargin == 1
    varargout{1} = subsref(defaults_suit, subs);
else
    defaults_suit = subsasgn(defaults_suit, subs, varargin{1});
end
