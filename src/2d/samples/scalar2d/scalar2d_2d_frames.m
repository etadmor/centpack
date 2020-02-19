function scalar2d_2d_frames(J,K,N)
%scalar2d_2d_frames   Generates snapshots of solution of the 2D scalar conservation law 
%   
%   scalar2d_2d_frames(J,K,N) where J and K are the number of grid cells along 
%   the x- and y-direction and N is the total number of outputs generated 
%   by CentPack, loads the output files of CentPack's scalar2d_2d_SD3 example 
%   and plots countours of the solution at dt_out intervals over the length 
%   of the simulation
%	
%	CentPack's output is written to the directory
%
%	CP_root/samples/scalar2d_2d_FD2/u_files/
%
%	where CP_root stands for your CentPack installation directory.  The 
%   data is loaded into matlab with the built-in load command and ploted 
%	over the solution domain, the resulting images are written as a .png 
%	files to the directories
% 
%	CP_root/samples/euler_2d_FD2/u_frames/
%	
%	A number of open source tools is available to create an animated 
%   sequence of the generated frames, one possiblitiy in UNIX-like systems 
%   is to run the following commands from the command window
%
%	> cd CP_root/samples/euler_2d_FD2/u_frames/
%	> convert -adjoin -delay 5 *.png u_movie.gif
%
%	These will generate the animation u_movie.gif.
%
%	Remark: convert is a command line application of the poen source graphics
%	suite ImageMagick(C), commonly distributed with UNIX-like systems
%	
%	Copyright 2004-2010 Jorge Balbas and Eitan Tadmor 
%   $Revision: 1.0 $  $Date: 2010/04/14
%

dx = 2*pi/J;
dy = 2*pi/K;

x = zeros(1,J);
y = zeros(1,K);

x(1) = 0.5*dx;
y(1) = 0.5*dy;

for j = 2:J
	x(j) = x(j-1) + dx;
end;

for k = 2:K
	y(k) = y(k-1) + dy;
end;

for n = 0:N-1

	count = int2str(n);
	s_u = strcat('u_files/u_', count);
	
	if n < 10
		count = strcat('0', count);
	end;
	
	S_u = strcat('u_frames/u_', count);
	
	u = load(s_u);
	
	contour(x,y,u',20);
    axis('square');
	print ('-dpng', '-r0', S_u);

end;
