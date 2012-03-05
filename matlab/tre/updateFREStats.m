function [fremean frecov frewin] = updateFREStats(newfre, fremean, frecov, frewin)
%% Written by Andrew D. Wiles, 2009-12-03.
%%
%% Here we update the FRE statistics for a single marker.
%% Input:
%%      newfre: 1x3 vector with the FRE components.
%%      fremean: 1x3 vector with the mean computed on the last frame.
%%      frecov: 3x3 covariance matrix computed on the last frame.
%%      frewin: the sliding window that holds the N FRE measurements that
%%      estimate the statistics.
%%  
%%  Here we add the newest data and remove the oldest.  Two methods are
%%  coded.  The first method assumes that the FRE mean is zero and hence
%%  the update can be reduced to adding the latest and removing the
%%  oldest. The second method assumes that the mean is not zero and that we
%%  need to recompute the covariance using the mean.  In simulation this
%%  works quite well but it still needs to be proven effective for real
%%  data.

winsize = size(frewin,1);

% hard code the method for the covariance computation.
covmethod = 1; % use the assumption that the FRE mean is zero.
%covmethod = 2; % use the matlab function cov.

switch covmethod
    case 1
        %update the mean.
        fremean = fremean + 1/winsize* (newfre - frewin(end,:));
        %update the covariance
        frecov = frecov + 1/(winsize-1) * (newfre'*newfre - frewin(end,:)'*frewin(end,:));
        %update the window.
        frewin = circshift(frewin,1);
        frewin(1,:) = newfre;
    case 2
        %update the window.
        frewin = circshift(frewin,1);
        frewin(1,:) = newfre;
        fremean = mean(frewin);
        frecov = cov(frewin); 
    otherwise
        error('Invalid covariance method given.');
end