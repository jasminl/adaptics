function varargout = roc(varargin)
% VL_ROC Compute the ROC curve
%  [TP,TN] = VL_ROC(Y, SCORE) computes the VL_ROC curve of the specified
%  data. Y are the ground thruth labels (+1 or -1) and SCORE is the
%  discriminant score associated to the data by a classifier (higher
%  scores correspond to positive labels).
%
%  [TP,TN] are the true positive and true negative rates for
%  incereasing values of the decision threshold.
%
%  [TP,TN,INFO] = VL_ROC(...) returns the following additional
%  informations:
%
%  INFO.EER::       Equal error rate.
%  INFO.AUC::       Area under the VL_ROC (AUC).
%  INFO.UR::        Uniform prior best op point rate.
%  INFO.UT::        Uniform prior best op point threhsold. 
%  INFO.NR::        Natural prior best op point rate.
%  INFO.NT::        Natural prior best op point threshold.
%
%  VL_ROC(...) plots the VL_ROC diagram in the current axis.
%
%  About the VL_ROC curve::
%    Consider a classifier that predicts as positive al lables whose
%    SCORE is not smaller than a threshold. The VL_ROC curve represents
%    the performance of such classifier as the threshold r is
%    varied. Denote:
%
%      P  = num of positive samples
%      N  = num of negative samples
%      TP = num of samples that are correctly classified as positive
%      TN = num of samples that are correctly classified as negative
%      FP = num of samples that are incorrectly classified as positive
%      FN = num of samples that are incorrectly classified as negative
%
%    Consider also the rates:
%
%                TP_ = TP / P,      FN_ = FN / P,
%                TN_ = TN / N,      FP_ = FP / N.
%
%    Notice that:
%
%                 P = TP  + FN ,    N = TN  + FP,
%                 1 = TP_ + FN_,    1 = TN_ + FP_.
%
%    The VL_ROC curve is the parametric curve (TP_, TN_) obtained
%    as the classifier threshold is changed.
%
%    The VL_ROC curve is contained in the square with vertices (0,0)
%    The (average) VL_ROC curve of a random classifier is a line which
%    connects (1,0) and (0,1).
% 
%    The VL_ROC curve is independent of the prior probability of positive
%    PPOS and negative labels PNEG. For instance, the empirical
%    expected error (01-risk) is
%
%        ERR = FP_ PPOS + FN_ PNEG,   PPOS = P/(P+N), 
%                                     PNEG = N/(P+N). 
%
%    An OPERATING POINT is a point on the VL_ROC curve, corresponding to
%    a certain threshold. Each operating point minimizes the empirical
%    error for certain label priors PPOS and PPNEG.  VL_ROC() computes
%    the following operating points:
%
%     Natural operating point:: Assumes PPOS = P/(P+N).
%     Uniform operating point:: Assumes PPOS = 1/2.
%
%  See also:: VL_HELP().
[varargout{1:nargout}] = vl_roc(varargin{:});
