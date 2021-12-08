%% first we're going to find the singularity at about 40 au (barycentric)
do_this = 1;
if do_this == 1
    clc
    amin = 40;
    amax = 41;
    na = 1000;
    na = na + 1;
    planetdir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3';
    planetfile = 'sun_and_planets_Horizons_barycentric_20220122.csv';
    savedir = planetdir;
    codedir = planetdir;
    yr = '20220122';
    A = laplace_plane_filesaver(amin,amax,na,planetdir,planetfile,savedir,codedir,yr);
    cd(codedir)
    a_vec = A(:,1);
    p0_vec = A(:,2);
    [~,ind1] = min(p0_vec);
    [~,ind2] = max(p0_vec);
    center_ind = floor( (ind1+ind2)/2 );
    a_bary_40 = a_vec(center_ind);
    disp('done')
end
%% now we're going to find the singularity at about 35 au (barycentric)
do_this = 1;
if do_this == 1
    clc
    amin = 34;
    amax = 36;
    na = 1000;
    na = na + 1;
    planetdir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3';
    planetfile = 'sun_and_planets_Horizons_barycentric_20220122.csv';
    savedir = planetdir;
    codedir = planetdir;
    yr = '20220122';
    A = laplace_plane_filesaver(amin,amax,na,planetdir,planetfile,savedir,codedir,yr);
    cd(codedir)
    a_vec = A(:,1);
    p0_vec = A(:,2);
    [~,ind1] = min(p0_vec);
    [~,ind2] = max(p0_vec);
    center_ind = floor( (ind1+ind2)/2 );
    a_bary_35 = a_vec(center_ind);
    disp('done')
end
%% save bin boundaries to file
do_this = 1;
if do_this == 1
    clc
    amin_template = [0 0 42 43 44 45 50 50 30 0 0 42 45 42]; % blanks at 1,2,10,11
    amax_template = [0 42 43 44 45 50 80 150 150 50 150 45 48 48]; % blank at 1
    amin_bary = amin_template;
    amin_bary([1 2 10 11]) = [a_bary_35 a_bary_40 a_bary_35 a_bary_35];
    amax_bary = amax_template;
    amax_bary(1) = a_bary_40;
    A = [amin_bary; amax_bary];
    A = transpose(A);
    savedir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3';
    codedir = savedir;
    cd(savedir)
    filename = 'amin_amax_bary_20220122.txt';
    writematrix(A,filename);
    cd(codedir)
    disp('done')
end
%% load bin boundaries from file, save Laplace plane for all bins (barycentric)
do_this = 1;
if do_this == 1
    clc
    tic
    datetime
    binsdir = '/Users/iggymatheson/Documents_off_iCloud/mm22_v3';
    binsfile = 'amin_amax_bary_20220122.txt';
    savedir = binsdir;
    codedir = binsdir;
    planetfile = 'sun_and_planets_Horizons_barycentric_20220122.csv';
    cd(binsdir)
    A = load(binsfile);
    cd(codedir)
    amin_vec = A(:,1);
    amax_vec = A(:,2);
    nbins = length(amin_vec);
    for ibin = 1:nbins
%     for ibin = 14:14
        amin = amin_vec(ibin);
        amax = amax_vec(ibin);
        na = ceil( (round(amax,3)-round(amin,3))*500 );
        disp([ibin nbins na])
        yr = '20220122';
        A = laplace_plane_filesaver(amin,amax,na,binsdir,planetfile,savedir,codedir,yr);
        cd(codedir)
    end
    toc
    disp('done')
end