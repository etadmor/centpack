function euler_2d_frames(J,K,N)
%euler_2d_frames   Generates snapshots of solution of 2D Euler's EQNs 
%   
%   euler_2d_frames(J,K,N) where J and K are the number of grid cells along 
%   the x- and y-direction and N is the total number of outputs generated 
%   by CentPack, loads the output files of CentPack's euler_2d_FD2 example 
%   and plots the density and pressure countours, and the velocity field 
%   at dt_out intervals over the length of the simulation
%	
%	CentPack's output is written to the directory 
%
%	CP_root/samples/euler_2d_FD2/rho_files/
%   CP_root/samples/euler_2d_FD2/u1_files/
%   CP_root/samples/euler_2d_FD2/u2_files/
%   CP_root/samples/euler_2d_FD2/p_files/
%
%	where CP_root stands for your CentPack installation directory.  The 
%   data is loaded into matlab with the built-in load command and ploted 
%	over the solution domain, the resulting images are written as a .png 
%	files to the directories
% 
%	CP_root/samples/euler_2d_FD2/rho_frames/
%	CP_root/samples/euler_2d_FD2/u_frames/
%	CP_root/samples/euler_2d_FD2/p_frames/
%	
%	A number of open source tools is available to create an animated 
%   sequence of the generated frames, one possiblitiy in UNIX-like systems 
%   is to run the following commands from the command window
%
%	> cd CP_root/samples/euler_2d_FD2/VAR_frames/  (VAR = rho, u, p)
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

dx = 1.0/J;
dy = 1.0/K;

x(1) = 0.5*dx;
y(1) = 0.5*dy;

for j = 2:J
    x(j) = x(j-1) + dx;
end

for k = 2:K
    y(k) = y(k-1) + dy;
end

xx = x(1:5:J);
yy = y(1:5:K);

for n = 0:N-1

	count=int2str(n);
	s_rho=strcat('rho_files/rho_', count);
	s_u1=strcat('u1_files/u1_', count);
	s_u2=strcat('u2_files/u2_', count);
	s_p=strcat('p_files/p_', count);
	
	if n<10
		count = strcat('0', count);
	end;
	
	S_rho = strcat('rho_frames/rho_', count);
	S_u = strcat('u_frames/u_', count);
	S_p = strcat('p_frames/p_', count);
	
	rho = load(s_rho);
	u1 = load(s_u1);
	u2 = load(s_u2);
	p = load(s_p);
	
	ux = u1(1:5:J,1:5:K);
	uy = u2(1:5:J,1:5:K);
	
	contour(x,y,rho',20);
    axis('square')
	print ('-dpng', '-r0', S_rho);
	
	quiver(xx,yy,ux',uy',1.5);
	axis([0 1 0 1]);
    axis('square')
	print ('-dpng', '-r0', S_u);
	
	contour(x,y,p',20);
    axis('square')
	print ('-dpng', '-r0', S_p);

end;
