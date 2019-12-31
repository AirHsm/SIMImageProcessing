function illumination = generate_illumination(img, k,num_rotations, num_phases)
	m = size(img, 1);
	cent = fix(m/2+1);
	x = linspace(1,m,m);
	y = linspace(1,m,m);
	[X,Y] = meshgrid(x,y);
	illumination = zeros(m,m,num_rotations*num_phases);

	n = 1;
	rotation = pi/20;
	for t = 1 : num_rotations
		phase = 0;
		k_x = k * cos(rotation);
		k_y = k * sin(rotation);
		kv = [k_x,k_y];
		rotation = rotation + pi/num_rotations;

		for p = 1:num_phases
			%illumination(:,:,n) = 0.5.*(1.0+cos(2*pi/m.*(k_x.*(X-cent)+k_y.*(Y-cent))+phase));
			illumination(:,:,n) = get_pattern(kv,phase,m,0.8,1);
			phase = phase + 2*pi/num_phases;
			n = n+1;
		end
	end

end
