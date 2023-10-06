function [data,time]=ah_get_data(filename,irec,tt,clp1,clp2,chp1,chp2)
%AH_GET_DATA Function to extract data from an ahx file, and (optionally) to
%filter it.
%
%This function is intended to emulate TSTOOL in its operation. As a result,
%it only performs internal calculations using positive-frequency Fourier
%components; a Hermitian spectrum is then generated before transforming
%back to the time domain. This has an additional advantage, of
%approximately halving the computational time required to execute the
%function.
%
%[data,time]=AH_GET_DATA(filename,irec), where filename is the path to a
%valid .ahx file and irec is a record number, returns the raw data contained
%in the record. No filtering is done in this case. The array 'time' will
%contain the time in seconds of each sample.
%
%[data,time]=AH_GET_DATA(filename,irec,tt) returns processed data. The
%instrument response has been deconvolved, and the data is convolved with
%an SRO instrument response. A time-domain cosine taper of length tt
%seconds has been applied. Note that SRO-response data with no time-domain
%taper may be obtained if tt=0.
%The filtering in this case is equivalent to the TSTOOL command:
% 'fi -sro -tt tt'
%
%[data,time]=AH_GET_DATA(filename,irec,tt,clp1,clp2) returns data filtered
%with a cosine low-pass filter between clp1 seconds and clp2 seconds (i.e.
%no periods less than clp1 will be passed). Time-domain tapering and an SRO
%response are also applied.  This is equivalent to the TSTOOL comand:
% 'fi -sro -tt tt -clp clp1s clp2s'
%
%[data,time]=AH_GET_DATA(filename,irec,tt,clp1,clp2,chp1,chp2) returns data
%filtered as above, but with the addition of a cosine high-pass filter
%characterised by chp1 seconds and chp2 seconds. This is equivalent to:
% 'fi -sro -tt -clp clp1s clp2s -chp chp1s chp2s'


%Check that number of arguments is permissible
disp(nargchk(2,7,nargin)); % nargin=2: No filters
                           %        3: SRO;tapers
                           %        5: SRO;tapers;CLP
                           %        7: SRO;tapers;CLP;CHP


if nargin==4 | nargin==6
    error('ah_get_data:argswrong','Impermissible number of arguments provided')
end
    
                           
%Open file
fid=fopen(filename,'r');
if fid<0
    error('ah_get_data:filenotfound','File not found. Please check that filename is correct, and provide path to file if not in current Matlab working directory')
end

fseek(fid,0,-1);
AAHNDT=165;
LAHHDR=270;
AAHDSN=12;
AAHDEL=166;
nrecs=0;
ibyte=AAHNDT*4;


%Find record
while (nrecs<irec && fseek(fid,ibyte,0)==0)
    npoints=fread(fid,1,'uint32','ieee-be');
    nrecs=nrecs+1;
    ibyte=(LAHHDR-1+npoints)*4;
end

%Get response function, even if we  won't use it - saves messing around
%with conditional fseeks

fseek(fid,(AAHDSN-AAHNDT-1)*4,0);
rfn=fread(fid,122,'float','ieee-be');


%Get sampling interval
fseek(fid,(AAHDEL-AAHDSN-122)*4,0);
smpin=fread(fid,1,'float','ieee-be');


%Read data
fseek(fid,(LAHHDR-AAHDEL-1)*4,0);
data=fread(fid,npoints,'float','ieee-be');
time=0:smpin:smpin*(npoints-1);

%Remove average
%data=remavl(npoints,data); % TOT

%If we just want data, we're done; otherwise, need to FFT
if nargin>2
    %Time-domain tapering
    ntaperp=ceil(tt/smpin);
    if ntaperp>npoints/2
        ntaperp=floor(npoints/2);
    end
    
    for i=1:ntaperp
        fac=0.5*(1-cos((i-1)*pi/ntaperp));
        data(i)=data(i)*fac;
        data(npoints+1-i)=data(npoints+1-i)*fac;
    end

    %FFT
    nspec=2^(nextpow2(npoints));
    maxspecpt=(nspec/2)+1; %Index of highest positive-frequency component (Nyquist component)
    spec=fft(data,nspec);
    %Only work with positive frequencies...
    f=0:1/(nspec*smpin):1/(2*smpin);

    %Divide by response; multiply by SRO
    for i=1:maxspecpt
        spec(i)=spec(i)*5000*1e6*sroins(2*pi*f(i));
        spec(i)=spec(i)/resp(2*pi*f(i),rfn);
    end

    if nargin>4 %CLP
        fl1=1/clp1;
        fl2=1/clp2;
        if nargin>6 %CHP
          fh1=1/chp1;
          fh2=1/chp2;
        end
        for i=1:maxspecpt-1
            fac=1;
            if f(i)>fl2
                fac=fac*0;
            else
                if f(i)>fl1
                    fac=fac*0.5*(1-cos((f(i)-fl2)*pi/(fl1-fl2)));
                end
            end
            if nargin>6
                if f(i)<fh1
                    fac=fac*0;
                else
                    if f(i)<fh2
                        fac=fac*0.5*(1-cos((fh1-f(i))*pi/(fh1-fh2)));
                    end
                end
            end
            spec(i)=spec(i)*fac;
        end
    end

    
    %Make spectrum Hermitian
    j=maxspecpt-1;
    for i=maxspecpt+1:nspec
        spec(i)=conj(spec(j));
        j=j-1;
    end
    
    %IFFT
    spec(1)=0;
    spec((nspec/2)+1)=0;
    invfft=ifft(spec,nspec);
    for i=1:npoints
        data(i)=invfft(i);
    end
end


fclose(fid);
end

function romega=resp(omega,rfn)
    romega=rfn(1)*rfn(2);
    i=0;
    nscale=0;
    while i<rfn(5)
        romega=romega*complex(-rfn(9+(4*i)),omega-rfn(10+(4*i)));
        i=i+1;
        if abs(romega)>1e20
            romega=romega*1e-20;
            nscale=nscale+1;
        end
        if abs(romega)<1e-20
            romega=romega*1e20;
            nscale=nscale-1;
        end
    end
    i=0;
    while i<rfn(3)
        romega=romega/complex(-rfn(7+(4*i)),omega-rfn(8+(4*i)));
        i=i+1;
        if abs(romega)>1e20
            romega=romega*1e-20;
            nscale=nscale+1;
        end
        if abs(romega)<1e-20
            romega=romega*1e20;
            nscale=nscale-1;
        end
    end
    romega=romega*((1e20)^nscale);
end

function sro=sroins(omega)
    s=complex(0,omega);
    z=(s^5)*(s+complex(0.1256,0))*(s+complex(50.1,0))*(s+complex(0,1.053))*(s+complex(0,-1.053));
    p=(s+complex(4.648,3.463))*(s+complex(4.648,-3.463))*(s+complex(0.1179,0.0))*(s+complex(40.73,0.0))*(s+complex(100.,0.0))*...
        (s+complex(0.1500,0.0))*(s+complex(264.,0.0))*(s+complex(3.928,0.0))*(s+complex(0.2820,0.0))*(s+complex(0.2010,0.2410))*...
        (s+complex(0.2010,-0.2410))*(s+complex(0.1337,0.1001))*(s+complex(0.1337,-0.1001))*(s+complex(0.0251,0.0))*(s+complex(0.00924,0.0))...
        *(s+complex(0.3334,0.2358))*(s+complex(0.3334,-0.2358))*(s+complex(0.1378,0.5873))*(s+complex(0.1378,-0.5873));
    sro=complex(9.237E+03,0.0)*(z/p);
end

function rdata=remavl(n,data)
    sumx=n*(n+1)/2;
    sumx2=sumx*(2*n+1)/3;
    sumy=sum(data(1:n));
    sumxy=sum((1:n)'.*data(1:n));
    det=(n*sumx2)-(sumx*sumx);
    if det ~= 0
        a=(n*sumxy-sumx*sumy)/det;
        bb=(sumx2*sumy-sumxy*sumx)/det;
    else
        a=0;
        bb=0;
    end
    for i=1:n
        rdata(i)=data(i)-i*a-bb;
    end
    
        
        
    
end