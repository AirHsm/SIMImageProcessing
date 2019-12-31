function mi = H(structure, illumination, otf)
	mi = zeros(size(illumination));
	if length(size(illumination))==2
		mi = abs(ifft2(otf .* fftshift(fft2(structure.*illumination))));
	else
		for n = 1:size(illumination,3)
			si = structure .* illumination(:,:,n);
			SI = fftshift(fft2(si));
			mi(:,:,n) = abs(ifft2(SI .* otf));
		end
	end
end
