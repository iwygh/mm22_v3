function A = laplace_plane_filesaver(amin,amax,na,planetdir,planetfile,savedir,codedir,yr)
a_vec = linspace(amin,amax,na);
a_vec = transpose(a_vec);
na = length(a_vec);
p0_vec = zeros(na,1);
q0_vec = zeros(na,1);
i0_vec = zeros(na,1);
Omega0_vec = zeros(na,1);
for ia = 1:na
    a_sample = a_vec(ia);
    cd(codedir);
    [p0,q0,i0,Omega0] = laplace_plane_function(a_sample,planetdir,planetfile);
    p0_vec(ia) = p0;
    q0_vec(ia) = q0;
    i0_vec(ia) = i0;
    Omega0_vec(ia) = Omega0;
%     disp([ia na])
end
i0_vec = i0_vec*180/pi;
Omega0_vec = Omega0_vec * 180/pi;
a_vec = real(a_vec);
i0_vec = real(i0_vec);
Omega0_vec = real(Omega0_vec);
A = [a_vec p0_vec q0_vec i0_vec Omega0_vec];
prefix_string = 'a_p0q0_i0Omega0_';
filename = horzcat(prefix_string,num2str(amin),'_',num2str(amax),'_',yr,'.txt');
if not(isfolder(savedir))
    mkdir(savedir)
end
cd(savedir)
writematrix(A,filename)
end
