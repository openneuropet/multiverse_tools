function [W, p, q] = STAPLE(varargin)

% Vecorized MATLAB implementation of the STAPLE algorithm by Warfield et al.
% Andreas Husch % 2013-07-10, 2015-07-06 % mail@andreashusch.de
%
% This code is an adaptation of the original code to work directly with segmented image from SPM
% Function:  [W, p, q] = STAPLE(P)
%
% IMPUT : P is a cell array of segmented image names (transfomed in D, data matrix of segmentations, dimensions VOXELS x IMAGES)
%         threshold to apply to W to get an estimated ground truth segmentation (W >= .5 by default)
% OUTPUT: W, est. weight matrix for each voxel
%         p, est. vector of sensitivities for each expert
%         q, est. vector of specificities for each expert
%         staple_mask, a binary image writen on the drive
%
% Literature: Warfield, Simon K., Kelly H. Zou, and William M. Wells. 
%            "Simultaneous truth and performance level estimation (STAPLE): 
%            an algorithm for the validation of image segmentation." 
%           Medical Imaging, IEEE Transactions on 23.7 (2004): 903-921.
%
% edits to work with images etc by Cyril Pernet. requires SPM in the path

%% check inputs
if nargin == 0
    help STAPLE
    return
else 
    P = varargin{1};
    if ~iscell(P)
        error('argument in must be a cell array of names')
    end
    
    if nargin == 1
        threshold = 0.5;
    else
        threshold = varargin{2};
    end
    clear varargin
end

if nargin > 2
    warn('more arguments in than needed??? only using arguments 1 and 2')
end

%% deal with images
if size(P,1) == 1; P=P'; end
N = size(P,1);
dim = NaN(N,3);
for img = 1:N
    if iscell(P)
        V(img) = spm_vol(P{img});
    else
        V(img) = spm_vol(P(img,:));
    end
    dim(img,:) = V(img).dim;
end 
    
% if sum(sum(dim==dim(1,:)) == N) ~= 3
%     error('imput images do not all have the same dimensions');
% end

%% make the matrix D
M = prod(dim(1,:));
D = NaN(M,N);
for img = 1:N
    d = spm_read_vols(V(img));
    d = d>0; % just in case not binary
    D(:,img) = d(:);
end
D = double(D);    

%% STAPLE-Algorithm by Warfield et al. for binary segementations
% --> Andreas Husch's code
     
    %% Parameters
    MAX_ITERATIONS = 30;
    EPSILON = 0.00001; % convergence criterion
    
    % Initial sensitivity and specificity parameter p(j),q(j) for all
    % raters j
    p(1:N) = 0.99999; 
    q(1:N) = 0.99999;
    Tprior = (sum(D(:))/length(D(:))); % note dependence on (sub)volume size, final result depends on this prior (which is not an implementation issue but a basic limitation of the EM approach)

    avgW = 1;
    W = zeros(1,length(D));
    
    %% EM
    
    for step=1:MAX_ITERATIONS 
        
        % E-Step
        % The following vectorized code is equivalent to this loop by MUCH
        % faster on current CPUs
        %     for i = 1:length(D)
        %         W(i) = ((prod(p(D(i,:))) * prod(1 - p(~D(i,:)))) * Tprior) / ... 
        %                ((prod(p(D(i,:))) * prod(1 - p(~D(i,:)))) * Tprior + (prod(q(~D(i,:))) * prod(1 - q(D(i,:))))) * (1- Tprior) ;
        %         %NOTE that prod([]) = 1
        %     end
    
        
        P = repmat(p,length(D), 1);
        Q = repmat(q,length(D), 1);
        P_given_D = P .* D; %TODO: use bsxfun instead of repmat?
        P_given_D(P_given_D(:)== 0) = 1; %
        Q_given_D = 1 - Q .* D;
        Q_given_D(Q_given_D(:)== 0) = 1; % alternative: initialize with 1 and set Q_given_D(D) = 1- P(D) 
        compP_given_not_D  = 1 - P .* ~D;
        compP_given_not_D(compP_given_not_D(:)== 0) = 1;
        compQ_given_not_D  = Q .* ~D;
        compQ_given_not_D(compQ_given_not_D(:)== 0) = 1;

        % W(i) can be interpretated as the prob. of voxel i beeing true (i.e. is part of the groundtruth y) for given p(1:N), q(1:N) 
        W = (prod(P_given_D') .* prod(compP_given_not_D') * Tprior) ./ ...
           ((prod(P_given_D') .* prod(compP_given_not_D') * Tprior) + (prod(Q_given_D') .* prod(compQ_given_not_D') * (1 - Tprior)));

        % Convergence?
        if(abs(avgW - sum(W) / length(W)) < EPSILON)
            break;
        end
        avgW = sum(W) / length(W);

        % M-Step
        p = (W * D) / sum(W(:)); % W * D = sum(W(D))
        q = ((1 -  W) * ~D) / sum(1 - W(:)); 
    end

%% write the result image
R = V(1);
location = fileparts(V(1).fname);
R.fname = [location filesep 'staple_mask.nii'];
R.descrip = 'STAPLE from binary images';
R.private.descrip = R.descrip;
Img = reshape(W,dim(1,:));
spm_write_vol(R,Img);
Img = Img.*(Img>threshold);
R.fname = [location filesep 'tstaple_mask.nii'];
spm_write_vol(R,Img);
