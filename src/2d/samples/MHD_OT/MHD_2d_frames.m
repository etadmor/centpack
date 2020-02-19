function MHD_2d_frames(J,K,N)
%MHD_2d_frames   Generates snapshots of solution of 2D MHD EQNs 
%   
%   MHD_2d_frames(J,K,N) where J and K are the number of grid cells along 
%   the x- and y-direction and N is the total number of outputs generated 
%   by CentPack, loads the output files of CentPack's MHD_2d_SD2 example 
%   and plots the density and pressure countours, and the velocity and 
%   magnetic vector fields at dt_out intervals over the length of the 
%   simulation
%	
%	CentPack's output is written to the directory 
%
%	CP_root/samples/MHD_2d_SD2/rho_files/
%   CP_root/samples/MHD_2d_SD2/u1_files/
%   CP_root/samples/MHD_2d_SD2/u2_files/
%   CP_root/samples/MHD_2d_SD2/b1_files/
%   CP_root/samples/MHD_2d_SD2/b2_files/
%   CP_root/samples/MHD_2d_SD2/p_files/
%
%	where CP_root stands for your CentPack installation directory.  The 
%   data is loaded into matlab with the built-in load command and ploted 
%	over the solution domain, the resulting images are written as a .png 
%	files to the directories
% 
%	CP_root/samples/MHD_2d_SD2/rho_frames/
%	CP_root/samples/MHD_2d_SD2/u_frames/
%	CP_root/samples/MHD_2d_SD2/b_frames/
%	CP_root/samples/MHD_2d_SD2/p_frames/
%	
%	A number of open source tools is available to create an animated 
%   sequence of the generated frames, one possiblitiy in UNIX-like systems 
%   is to run the following commands from the command window
%
%	> cd CP_root/samples/euler_2d_FD2/VAR_frames/  (VAR = rho, u, b, p)
%	> convert -adjoin -delay 5 *.png VAR_movie.gif
%
%	These will generate the animation VAR_movie.gif.
%
%	Remark: convert is a command line application of the poen source graphics
%	suite ImageMagick(C), commonly distributed with UNIX-like systems
%	
%	Copyright 2004-2010 Jorge Balbas and Eitan Tadmor 
%   $Revision: 1.0 $  $Date: 2010/04/14
%

x = zeros(1,J);
y = zeros(1,K);

dx = 2*pi/J;
dy = 2*pi/K;

x(1) = 0.5*dx;
y(1) = 0.5*dy;

for j = 2:J
    x(j) = x(j-1) + dx;
end

for k = 2:K
    y(k) = y(k-1) + dy;
end

xx = x(1:4:J);
yy = y(1:4:K);

for n = 0:N-1

% The following assumes that the numerical results for density, velocity field
% components, etc, that your code generates are in directories called rho_files,
% u1_files, u2_files, etc., with names rho_0, rho_1, rho_23, ..., rho_60 (note
% the single digit in rho_0, rho_1, graphic utilities need that changed to 
% rho_00, ..., rho_09, that's taken care below

	% converts the counter to a character
	
    count=int2str(n);

    % to read, recursively, the files in rho_files/, u1_files/, etc

    s_rho=strcat('rho_files/rho_', count);
    s_u1=strcat('u1_files/u1_', count);
    s_u3=strcat('u3_files/u3_', count);
    s_b1=strcat('b1_files/b1_', count);
    s_b3=strcat('b3_files/b3_', count);
    s_p=strcat('p_files/p_', count);

    % This fixes the rho_0 -> rho_00 conversion

    if n < 10
        count = strcat('0', count);
    end;

    % This will hold the names of the image files to be written.  You need to 
    % create the directories rho_movie, u1_movie, etc. in the same parent
    % directory where the rho_files/, u1_files/ and others are

    S_rho = strcat('rho_frames/rho_', count);
    S_u = strcat('u_frames/u_', count);
    S_b = strcat('b_frames/b_', count);
    S_p = strcat('p_frames/p_', count);

    % This automatically loads the numerical results one by one

    rho = load(s_rho);
    u1 = load(s_u1);
    u3 = load(s_u3);
    b1 = load(s_b1);
    b3 = load(s_b3);
    p = load(s_p);

    % This are coarser samples of the velocity and magnetic field components
    % to plot their vector fields

    ux = u1(1:4:J,1:4:K);
    uy = u3(1:4:J,1:4:K);

    bx = b1(1:4:J,1:4:K);
    by = b3(1:4:J,1:4:K);

    % This creates a color-filled contour plot of the density, and saves it as
    % rho_movie/rho_00

    if n > 0
        contourf(x,y,rho',20);
        axis('square');
        print ('-dpng', '-r0', S_rho);
    end

    quiver(xx,yy,ux',uy',1.5);
    axis([0 2*pi 0 2*pi]);
    axis('square');
    print ('-dpng', '-r0', S_u);
    %hold off;

    quiver(xx,yy,bx',by',1.5);
    axis([0 2*pi 0 2*pi]);
    axis('square');
    print ('-dpng', '-r0', S_b);
    %hold on;

    contourf(x,y,p',20);
    axis('square');
    print ('-dpng', '-r0', S_p);

end;
