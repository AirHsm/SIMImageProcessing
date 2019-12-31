% Fringe modulation factor estimation using In & Out focus information

function MfList = MfEstimate2(SI, WF, OTF, K, Phase)
  % Mf = MfEstimate2(si, wf)
  % addpath('../jRL_Deconvolution');
  [M,N,H] = size(SI);

  OTF = OTF / max(OTF(:));
  MfList = zeros(2,H);

  WF = jRL_deconvolve(WF, ones(M,N), OTF, 20);

  for i = 1:H
    si = SI(:,:,i);
    k = K(i,:);
    si = jRL_deconvolve(si, ones(M,N), OTF, 20);
    phi = Phase(i);
    MfOpt0 = @(a)MfOpt2(a,  WF, k, phi, si);
    options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
    Mf = fminunc(MfOpt0, 0.5, options);

    MfList(2,i) = Mf;
    MfList(1,i) = 1;
  end
