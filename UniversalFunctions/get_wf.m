function wf = get_wf(data)
    if length(size(data)) == 2
        wf = data;
    else
        m = size(data,1);
        wf = zeros(m,m);
        for i = 1:size(data,3)
            wf = wf+data(:,:,i);
        end
    end
    wf = img_normalize(wf);
    