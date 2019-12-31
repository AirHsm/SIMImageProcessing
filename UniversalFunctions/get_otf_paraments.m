function k = get_otf_paraments(otf)
	m = size(otf,1);
	cent = fix(m/2+1);
	l = otf(cent,:);
	for n =1:length(l)
		if l(n) > max(l)*1e-2
			break;
		end
	end
	k = cent-n;
end
