%%
function [r, z, t] = bendingfree_opt(r1, r2, N, NTF, ITF)
  %
  mu0 = 4 * pi * 1e-7;
  %
  tol = 1e-3;
  %
  % Parametrization
  %
  t = linspace(-pi/2, pi/2, N)';
  r = (r2 + r1)/2+(r2 - r1)/2 * cos(t + pi/2);
  %z = zeros(N, 1);
  z = (r2 - r1)/2 * sin(t + pi/2);
  %
  % Initial guess for z
  %
  ra = r;
  za = z;
  %
  errorr = 1;
  errorz = 1;
  iter = 0;
  %
  while errorz > tol || errorr > tol
    %
    iter = iter + 1;
    %
    % Generate input for ANSYS from inboard (r=r1) to outboard leg (r=r2)
    % Additional row to complete straight leg
    %
	rs = [r; r(N)];
    zs = [z; 0];
    %
    rs = [rs(end:-1:1); rs(1:end-1)];
    zs = [zs(end:-1:1); -zs(1:end-1)];
    figure(1)
    plot(rs,zs,'.');hold on
    axis equal
    ylabel('z [m]','fontsize',20,'Interpreter','latex')
    xlabel('r [m]','fontsize',20,'Interpreter','latex')
    drawnow
    %
    % dlmwrite('radialMatlab.txt', rs);
    % dlmwrite('verticalMatlab.txt', zs);
    fileID = fopen('radialMatlab_VNS.txt','w');
    fprintf(fileID, '%f \n',rs);
    fclose(fileID);
    fileID = fopen('verticalMatlab_VNS.txt','w');
    fprintf(fileID, '%f \n',zs);
    fclose(fileID);
    %
    % Compute field in Center Line / Call ANSYS
    fileID = fopen('parametriMatlab.txt','w');
    fprintf(fileID, '!!! Parametri da Matlab \n\n');
    fprintf(fileID, '\n n_TF =  %f \n I_TF = %f \n N_NODES =  %f',NTF,ITF,N);
    fclose(fileID);
    ! "C:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\ANSYS221.exe" -b -np 6 -i MAIN.dat -o MAIN.out -j MAIN
    delete(sprintf('main.lock'));     
    for i=0:5
        delete(sprintf('MAIN%G.err',i));     
        delete(sprintf('MAIN%G.out',i));     
    end
    % Btor = zeros(N+1, 1);
    % for i = 1:N+1
    %    Btor(i) = mfrz(rs(i), zs(i)); 
    % end
    %
    A = readmatrix('B_.csv');
    Btor = A(:,4).*2;
    Btor = Btor(2:N+1);
    rr = r(end:-1:1); 
    fig2 = figure(2);
    plot(rr,Btor);hold on 
    axis equal
    grid on
    drawnow

    Btor = 2 * pi * Btor / (mu0 * NTF * ITF);
    %
    % Integrate field from r1 to r, int(B) from r1 to r is converted into a
    % function to save computing time
    %
    intB = zeros(N, 1);
       
    %zr = z(end:-1:1);
    for i = 2 : N
      %intB(i) = intB(i-1) + 0.5 * (mfrz(rr(i-1), zr(i-1)) + mfrz(rr(i), zr(i))) * (rr(i) - rr(i-1));
      intB(i) = intB(i-1) + 0.5 * (Btor(i-1) + Btor(i)) * (rr(i) - rr(i-1));
    end
    intB_fh = @(ri) interp1(rr, intB, ri, 'linear', 'extrap');
    %
    % Compute tensile force in coil
    %
    T = mu0 * NTF * ITF^2 / 8 / pi * intB(end);
    k = 4 * pi * T / (mu0 * NTF * ITF^2);
    %
    % For each theta compute r(theta) using Newton method
    %
    x0 = r2;
    for i = 1 : N
      xx = x0;
      error = 1.0;
      while error > tol
        F = intB_fh(x0) + k * (sin(t(i)) - 1);
        % dF = mfrz(x0, z(i));
        % d2F = dmfrz_dr(x0, z(i)); % + dmfrz_dz(x0, z(i)) * tan(t(i));
        % xx = x0 - 2 * F * dF / (2 * dF^2 - d2F * F);
        xx = x0 - F / interp1(rr, Btor, x0, 'linear', 'extrap');
        error = abs(xx - x0) / abs(x0);
        x0 = xx;
      end
      r(i) = xx;
    end
    %
    % For each theta compute z(theta)
    %
    Btor = Btor(end:-1:1);
    for i = 2 : N
      % z(i) = z(i-1) - k * (sin(t(i-1)) / mfrz(r(i-1), z(i-1)) + sin(t(i)) / mfrz(r(i), z(i))) /2 * (t(i) - t(i-1));
      z(i) = z(i-1) - k * (sin(t(i-1)) / Btor(i-1) + sin(t(i)) / Btor(i)) /2 * (t(i) - t(i-1));
    end
    %
    % Compare new z(theta) with previous z(theta)
    %
    errorr = norm(r - ra) / norm(r);
    errorz = norm(z - za) / norm(z);
    %
    % Update previous z(theta)
    %
    ra = r;
    za = z;
    %
    fprintf('- Iter. %i - T = %f - r error = %f - z error = %f\n', iter, T*1e-6, errorr, errorz);
    %
  end
  % Complete with inboard straight leg
  %
  t = [t; pi/2] * 180 / pi;
  r = [r; r(N)];
  z = [z; 0];
  %
  t = t(end:-1:1);
  r = r(end:-1:1);
  z = z(end:-1:1);
  %
  rho = zeros(N+1, 1);
  Tc = zeros(N+1, 1);
  dzdr = zeros(N+1, 1);
  dzdr2 = zeros(N+1, 1);
  %
  Btor = Btor(end:-1:1);
  %
  for i = 3:N
    h = r(i+1) - r(i);
    t = r(i) - r(i-1);
    dzdr(i) = (z(i+1) - z(i-1)) / (h + t);
    dzdr2(i) = 2 * (z(i+1) + z(i-1) - 2 * z(i) - dzdr(i) * (h - t)) / (h^2 + t^2);
     
    rho(i) = (1 + dzdr(i)^2)^(3/2) / dzdr2(i);
    Tc(i) = -rho(i) * ITF / 2 * (Btor(i) * NTF * ITF * mu0 / 2 / pi);
  end
  %
  fig = figure(100);
  plot(r(3:end-1), Tc(3:end-1) * 1e-6, 'r', 'LineWidth', 2);
  % print(fig,'-dpng','-r300');
  % print(fig2,'-dpng','-r300');
  %
end
