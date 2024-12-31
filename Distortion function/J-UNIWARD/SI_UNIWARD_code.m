function SI_UNIWARD_code(cover_img,stego_img,payload,basecost)
%           EXAMPLE - USING "SI-UNIWARD" embedding distortion
%
% -------------------------------------------------------------------------
% Copyright (c) 2013 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Author: Vojtech Holub
% -------------------------------------------------------------------------
% Contact: vojtech_holub@yahoo.com
%          fridrich@binghamton.edu
%          http://dde.binghamton.edu
% -------------------------------------------------------------------------
%clc; clear all; close all;

% load precover
precover_path = strrep(cover_img,'jpg','pgm');
precover = imread(precover_path);

% set payload
payload = single(payload);
% set quality factor
quality_factor = int32(75);

%% Run default embedding
% distortion = SI_UNIWARD(precover, cover_path, stego_path, payload, quality_factor);

%% Embedding extraction can be configured in following way:
%
% Any number of the following settings can be included
%
%   field               type    default value  description
%   -----------------------------------------------------------------------
%   config.STC_h        uint32      0  	       0 for optimal coding simulator otherwise sets STC submatrix height (try 7-12)
%	config.seed         int32       0 	       random seed
%
%   Example:
%
config.STC_h = uint32(0);
% config.seed  = int32(123);

distortion = SI_UNIWARD(precover, cover_img, stego_img, payload, quality_factor, config);
       
end