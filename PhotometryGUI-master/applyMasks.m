function [ avgs ] = applyMasks(masks, img )
%APPLYMASKS Extract mask-wise averages from an image
%   masks   an MxNxD matrix where D is the number of masks
%   img     the MxN image to apply the masks to
[M, N, D] = size(masks);
avgs = zeros(1, D);
if D == 1
    avgs = mean(mean(img(masks)));
else
    for i = 1:D
        avgs(i) = mean(img(masks(:,:,i)));
    end
end

% %% Plot masks to check: Sometimes size(masks) ~= size(img)! Need to fix bug in ROI calculation
% figure;
% for i_m = 1:size(masks, 3)
%     subplot(3, 2, i_m);
%     imagesc(masks(:, :, i_m));
%     axis equal tight
% end
