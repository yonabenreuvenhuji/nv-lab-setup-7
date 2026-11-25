function ok = testFitPlots(stubDir)
%TESTFITPLOTS Run every plotting path of all three robust fits. ~40 s, no display needed.
%
%   testFitPlots                 % uses the real graphics
%   testFitPlots('/path/stubs')  % headless: shadow the graphics calls
%
% ===================================================================
% WHY A TEST FOR A PLOT
% ===================================================================
%
% The numbers are graded by testRamseyFit and friends. The DRAWING is not, and it is where
% the mistakes actually landed. In one editing session the plotting code acquired three
% separate crash bugs that no numerical test could see:
%
%   * two helper functions deleted while a duplicate was being removed, so the title block
%     called functions that no longer existed
%   * rabiFitRobust calling describeExp, which was never defined anywhere -- sitting in a
%     catch block, so it only fired on the day something else had already gone wrong
%   * o.periodogram indexed as a two-column matrix when it is a struct with .f and .sse
%
% Every one of those is an "undefined function" or "out of bound" at the moment you most
% want a figure. The last one is the clearest argument for this file: showScan is OFF by
% default, so that bug would have waited until someone turned it on to debug something else,
% and then broken in their hands.
%
% This walks every panel combination -- 42 of them across all five robust fits -- and only
% checks that nothing throws.
% It does not check that the figure looks right; eyes do that. It checks that there IS a
% figure.
%
% HEADLESS USE. Pass a directory holding no-op shadows of figure, axes, subplot, plot,
% errorbar, fill, area, title, xlabel, ylabel, legend, grid, box, hold, set, get, xlim,
% ylim, uistack, print, gca and gcf, and it runs with no display and no renderer. get must
% answer 'figure' for handle 1 and 'axes' otherwise, and xlim/ylim must return [0 1].

    if exist('stubDir', 'var') && ~isempty(stubDir); addpath(stubDir, '-begin'); end
    fprintf('\n=== testFitPlots ===\n');

    ok = true;
    % ESR, with and without the residual panel
    f=2810:1.5:2930; n=numel(f);
    y=1-0.019./(1+((f-2855)/8).^2)-0.019./(1+((f-2888)/7).^2)+3e-4*nz(n,1);
    for sr = [false true]
      try
        esrFitRobust(struct('f',f,'y',y,'e',3e-4*ones(1,n),'plotIt',true,'verbose',false,'showResid',sr));
        fprintf('  PASS  esr plot, showResid=%d\n', sr);
      catch e; fprintf('  FAIL  esr plot showResid=%d: %s\n', sr, e.message); ok=false; end
    end
    % RABI, with and without resid and scan
    t=0.02:0.01:1.0; n=numel(t);
    yr=1-0.06*exp(-t/1.2).*(1-cos(2*pi*13.5*t))/2 + 4e-4*nz(n,2);
    for sr = [false true]
      try
        rabiFitRobust(struct('t',t,'y',yr,'e',4e-4*ones(1,n),'plotIt',true,'verbose',false,'showResid',sr));
        fprintf('  PASS  rabi plot, showResid=%d\n', sr);
      catch e; fprintf('  FAIL  rabi plot showResid=%d: %s\n', sr, e.message); ok=false; end
    end
    try
      rabiFitRobust(struct('t',t,'y',yr,'e',4e-4*ones(1,n),'plotIt',true,'verbose',false,'showScan',true));
      fprintf('  PASS  rabi plot, showScan\n');
    catch e; fprintf('  FAIL  rabi showScan: %s\n', e.message); ok=false; end
    % RAMSEY: 1 and 2 projections x all four panel combinations
    A=2.158; dt=1/(A*15); t=0.005:dt:2.5; n=numel(t); amp=0.0032*0.9946/2;
    G=(cos(2*pi*(0.3-2.19)*t)+cos(2*pi*0.3*t)+cos(2*pi*(0.3+2.19)*t))/3; D=exp(-(t/1.22).^1.6);
    Y=[0.9946+amp*D.*G+1.3e-4*nz(n,3); 0.9948-1.15*amp*D.*G+1.3e-4*nz(n,4)];
    for np = [2 1]
      for sr = [false true]
        for sp = [false true]
          YY = Y(1:np,:); EE = 1.3e-4*ones(np,n);
          try
            ramseyFitRobust(struct('t',t,'y',YY,'e',EE,'plotIt',true,'verbose',false, ...
                                   'showResid',sr,'showSpectrum',sp));
            fprintf('  PASS  ramsey plot, %d proj, resid=%d spec=%d\n', np, sr, sp);
          catch e; fprintf('  FAIL  ramsey %dproj resid=%d spec=%d: %s\n', np,sr,sp,e.message); ok=false; end
        end
      end
    end
    % ECHO: monotone and revival data x 1 and 2 projections x the residual panel
    tc=linspace(1,90,100); n=numel(tc);
    for wr = [false true]
      D=exp(-(tc/45).^2.2);
      if wr; D=D.*ecomb(tc,18.7,0.15); end
      Ye=[0.99-0.013*D+1.5e-4*nz(n,7); 0.991+0.013*D+1.5e-4*nz(n,8)];
      Ye(:,20)=0.99;                    % a diluted point: exercises the flag and the broken line
      for np = [2 1]
        for sr = [false true]
          try
            echoFitRobust(struct('t',tc,'y',Ye(1:np,:),'e',1.5e-4*ones(np,n), ...
                                 'plotIt',true,'verbose',false,'showResid',sr));
            fprintf('  PASS  echo plot, revivals=%d %dproj resid=%d\n', wr, np, sr);
          catch e; fprintf('  FAIL  echo revivals=%d %dproj resid=%d: %s\n', wr,np,sr,e.message); ok=false; end
        end
      end
    end
    % the no-contrast path, where the envelope must be SUPPRESSED rather than drawn from an
    % amplitude that means nothing. It is the branch that only fires on a bad run, which is
    % exactly when nobody wants to discover it is broken.
    try
      echoFitRobust(struct('t',tc,'y',0.99+1.5e-4*nz(n,9),'e',1.5e-4*ones(1,n), ...
                           'plotIt',true,'verbose',false));
      fprintf('  PASS  echo plot, no measurable contrast (envelope suppressed)\n');
    catch e; fprintf('  FAIL  echo no-contrast plot: %s\n', e.message); ok=false; end

    % T1: stretched and two-component x 1 and 2 projections x the residual panel x linear/log x
    tt = unique(round(logspace(log10(20), log10(20000), 24))); nn = numel(tt);
    for dbl = [false true]
      if dbl; Dt = 0.4*exp(-tt/300) + 0.6*exp(-tt/4000); else; Dt = exp(-(tt/2500).^0.8); end
      Yt = [1 + 3e-4*nz(nn,11); 1 - 0.03*Dt + 3e-4*nz(nn,12)];
      Yt(2, 5) = 1.0;                    % a diluted point: the flag and the broken line
      md = 'stretched'; if dbl; md = 'auto'; end
      for np = [2 1]
        for sr = [false true]
          for lg = [true false]
            try
              t1FitRobust(struct('t',tt,'y',Yt(1:np,:),'e',3e-4*ones(np,nn),'plotIt',true, ...
                                 'verbose',false,'showResid',sr,'logX',lg,'model',md));
              fprintf('  PASS  t1 plot, double=%d %dproj resid=%d logX=%d\n', dbl, np, sr, lg);
            catch e
              fprintf('  FAIL  t1 plot double=%d %dproj resid=%d logX=%d: %s\n', dbl,np,sr,lg,e.message); ok=false;
            end
          end
        end
      end
    end
    % the no-decay path, where the plateau line and the T1 marker must be SUPPRESSED rather than
    % drawn from an amplitude that means nothing.
    try
      t1FitRobust(struct('t',tt,'y',1 + 3e-4*nz(nn,13),'e',3e-4*ones(1,nn), ...
                         'plotIt',true,'verbose',false));
      fprintf('  PASS  t1 plot, no measurable decay\n');
    catch e; fprintf('  FAIL  t1 no-decay plot: %s\n', e.message); ok=false; end

    % and the verbose reports, which is the other path my surgery touched
    try
      ramseyFitRobust(struct('t',t,'y',Y,'e',1.3e-4*ones(2,n),'plotIt',false,'verbose',true));
      fprintf('  PASS  ramsey verbose report\n');
    catch e; fprintf('  FAIL  ramsey report: %s\n', e.message); ok=false; end
    try
      D=exp(-(tc/45).^2.2).*ecomb(tc,18.7,0.15);
      echoFitRobust(struct('t',tc,'y',[0.99-0.013*D+1.5e-4*nz(n,7); 0.991+0.013*D+1.5e-4*nz(n,8)], ...
                           'e',1.5e-4*ones(2,n),'plotIt',false,'verbose',true,'T1us',400));
      fprintf('  PASS  echo verbose report (revival + T1 branches)\n');
    catch e; fprintf('  FAIL  echo report: %s\n', e.message); ok=false; end
    try
      Dt2 = exp(-(tt/2500).^0.8);
      t1FitRobust(struct('t',tt,'y',[1 + 3e-4*nz(nn,11); 1 - 0.03*Dt2 + 3e-4*nz(nn,12)], ...
                         'e',3e-4*ones(2,nn),'plotIt',false,'verbose',true,'model','auto'));
      fprintf('  PASS  t1 verbose report (control, reach and second-component branches)\n');
    catch e; fprintf('  FAIL  t1 report: %s\n', e.message); ok=false; end


    fprintf('\n');
    if ok
        fprintf('ALL PLOT PATHS RUN -- 42 panel combinations across five fitters, nothing thrown.\n\n');
    else
        fprintf('*** FAILURES ABOVE -- a figure will crash in your hands ***\n\n');
    end
end

function R = ecomb(t, Trev, wFrac)
% The revival comb, written out here rather than borrowed from echoFitRobust: a test that
% generates its data with the function under test cannot detect an error in that function.
    w = wFrac*Trev; K = floor(max(t)/Trev) + 2;
    R = zeros(size(t)); R0 = 0;
    for k = 0:K
        R  = R  + exp(-((t - k*Trev)/w).^2);
        R0 = R0 + exp(-((0 - k*Trev)/w).^2);
    end
    R = R / R0;
end

function r = nz(n, seed)
    m = 2^31 - 1;
    x = mod(1103515245*(seed + 12345) + 12345, m);
    u = zeros(1, 2*n);
    for i = 1:2*n
        x = mod(1103515245*x + 12345, m);
        u(i) = (x + 0.5)/m;
    end
    r = sqrt(-2*log(u(1:n))) .* cos(2*pi*u(n+1:2*n));
end
