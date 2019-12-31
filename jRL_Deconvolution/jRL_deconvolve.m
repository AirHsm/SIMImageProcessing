function result = jRL_deconvolve(sim_data, illumination, otf, num_iterations)
wf = get_wf(sim_data);
estimate = wf;
% hwait = waitbar(0,'jRL deconvolving...');
for iter = 1:num_iterations
	%expected_data = zeros(size(sim_data));
	expected_data = H(estimate, illumination, otf);
	ratio = sim_data ./ expected_data;
	estimate = estimate .* HT(ratio, illumination, otf) ./ HT(ones(size(illumination)), illumination, otf);
	estimate(isnan(estimate)) = 1;
    % str = ['jRL deconvolving...'];
    % waitbar(iter/num_iterations,hwait,str)

    % correct corners
    estimate(1,1) = mean([estimate(2,1), estimate(1,2), estimate(2,2)]);
    estimate(end,1) = mean([estimate(end-1,1), estimate(end,2), estimate(end-1,2)]);
    estimate(1,end) = mean([estimate(2,1), estimate(1,end-1), estimate(2,end-1)]);
    estimate(end,end) = mean([estimate(end-1,end), estimate(end,end-1), estimate(end-1,end-1)]);
end
% close(hwait);

result = estimate;
% result = result / max(result(:));
% fprintf('Done deconvolution!\n');
end
