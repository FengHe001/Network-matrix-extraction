function [Mask, ThresholdSelected] = f_IterativeThresholding(Im, MinThreshold, MaxThreshold, ThresholdStep, display)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
display = 0;
Steps = [MinThreshold : ThresholdStep : MaxThreshold];
for s = 1:size(Steps, 2)
   MaskThis = Im > Steps(s);
   AreaThis = sum(MaskThis(:));
   Areas(s) = AreaThis;
   if s > 1
        Deltas(s) = Areas(s) - Areas(s-1);
   end
end
if display
    figure
    plot(Areas)
end

ThresholdID = find(Deltas == min(Deltas));
ThresholdSelected = Steps(ThresholdID);
Mask = Im > ThresholdSelected;
% imtool(max(Mask,[],3))

end

