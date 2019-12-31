function structure = HT(mi, illumination, otf)
	structure = zeros(size(otf));
	if length(size(illumination)) == 2
		structure = abs(ifft2(otf.*fftshift(fft2(mi)))).*illumination;
	else
		for n = 1:size(illumination,3)
			MI = fftshift(fft2(mi(:,:,n)));
			structure = structure + abs(ifft2(MI .* otf).*illumination(:,:,n));
		end
	end
end
