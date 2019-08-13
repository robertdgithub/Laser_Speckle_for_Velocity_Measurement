%时空联合衬比分析方法（stLASCA）的衬比均值最大测量误差约5%
%而时间平均的空间衬比分析方法（sLASCA）和空间平均的时间衬比
%分析方法（tLASCA）的衬比均值最大测量误差为13%以上
% Objective:Let the newbie grasp how to use Laser speckle sk
% imaging as soon as possible
% Author:slni in IBP
% date:2005/07/01
% without permission, others have no right to modify the code!
% Speckle：显示原始散斑图
% SK：空间衬比图
% SV：由空间衬比得到的流速图
% TK：时间衬比图
% TV：由时间衬比得到的流速图
% stLASCA：时空联合衬比图
% sLASCA：时间平均的空间衬比图
% tLASCA：空间平均的时间衬比图
% nLASCA：标准化散斑衬比图
% addpt：考虑到静态散射成分的时间衬比图
% addpk：考虑到静态散射成分的空间衬比图
% AVS：动静脉分辨图

function varargout = SpecUI(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SpecUI_OpeningFcn, ...
    'gui_OutputFcn',  @SpecUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% --- Executes just before SpecUI is made visible.

function SpecUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.BlackThreshold,'value',0);
set(handles.WhiteThreshold,'value',4096);
handles.Cimage=0;
handles.bmroi=0;
handles.enableSK=0;
handles.enableSV=0;
handles.enableTK=0;
handles.enableTV=0;
handles.enableTSub=0;
handles.enablestLASCA=0;
handles.enabletLASCA=0;
handles.enablesLASCA=0;
handles.enablenLASCA=0;
handles.enableaddpt=0;
handles.enableaddpk=0;
handles.enableAVS=0;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes SpecUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function varargout = SpecUI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function Load_Callback(hObject, eventdata, handles)
SlideWindowDefault=5;
SlideWindow=str2num(get(handles.SlideWindow,'string'));
loadmodelindex=get(handles.LoadModel,'value');
framenum=str2double(get(handles.Frames,'string'));
groupnum=str2double(get(handles.Groups,'string'));
if(framenum>=1 && framenum<=1200 && groupnum>=1 && groupnum<=100)
    [fname pname]=uigetfile({'*.b16','16bit Orginal PixelFly Files(*.b16)';...
        '*.tif','tiff Data Files(*.tif)';...
        '*.tiff','tiff Data Files(*.tiff)';...
        '*.bmp','bmp Data Files(*.bmp)';...
        '*.jpg','jpeg Data Files(*.jpg)';...
        '*.b16;*.tif;*.bmp;*.jpg','Supportted file formats(*.b16,*.tif,*.tiff,*.bmp,*.jpg)'},...
        'Load Files for Processing');
    if fname==0
        errordlg('error: Fail to load the file name!','Load file problem');
    end
    handles.currentgroup=1;
    handles.FirstGroup=1;
    handles.LastGroup=groupnum;

    set(handles.Current,'Min',0,'Max',groupnum);
    handles.framenum=framenum;
    handles.groupnum=groupnum;

    handles.load_filetype=fname(strfind(fname,'.')+1:length(fname));

    intervIndex = strfind(fname,'_');
    if intervIndex >= 1
        intervIndex;
    else
        intervIndex=strfind(fname,'-');
    end
    intervChar = fname(intervIndex);

    filelen = length(fname(intervIndex+1:strfind(fname,'.')-1));

    firstgroupindex=str2num(fname(1:intervIndex-1));
    firstfileindex=str2num(fname(intervIndex+1:strfind(fname,'.')-1));

    if strcmp(handles.load_filetype,'b16')
        imagetemp=readb16(sprintf('%s%s',pname,fname));
    else
        imagetemp=imread(sprintf('%s%s',pname,fname));
    end
    handles.im=imagetemp;
    tmp=handles.im;
    tmp=double(tmp);
    a = spec2contr(tmp,SlideWindow);
    %load all the files to be processed
    [SpecX SpecY]=size(a);
    [OriFigX OriFigY]=size(imagetemp);

    SpeckleImg=zeros(SpecX,SpecY,framenum);
    handles.GroupSpecImg=zeros(OriFigX,OriFigY,groupnum);

    %     set(handles.checkSK,'Enable','off');
    if handles.enableSK==1
        handles.GroupKImg=zeros(SpecX,SpecY,groupnum);
    else
        set(handles.SK,'Enable','off');
    end
    
    %     set(handles.checkstLASCA,'Enable','off');
    if handles.enablestLASCA==1
        handles.GroupstVeloImg=zeros(OriFigX,OriFigY,groupnum);
    else
        set(handles.stLASCA,'enable','off');
    end
    
    %     set(handles.checksLASCA,'Enable','off');
    if handles.enablesLASCA==1
        handles.GroupsImg=zeros(SpecX,SpecY,groupnum);
    else
        set(handles.sLASCA,'enable','off');
    end
    
    %     set(handles.checknLASCA,'Enable','off');
    if handles.enablenLASCA==1
        handles.GroupnImg=zeros(SpecX,SpecY,groupnum);
    else
        set(handles.nLASCA,'enable','off');
    end
    
    %     set(handles.checktLASCA,'Enable','off');
    if handles.enabletLASCA==1
        handles.GrouptVeloImg=zeros(OriFigX,OriFigY,groupnum);
    else
        set(handles.tLASCA,'enable','off');
    end
    
    % 考虑到静态散射光
    if handles.enableaddpt==1
        handles.GroupaddVeloImg=zeros(OriFigX,OriFigY,groupnum);
    else
        set(handles.addpt,'enable','off');
    end
    
    % 考虑到静态散射时的空间衬比计算方法
    if handles.enableaddpk==1
        handles.GroupaddpkImg=zeros(SpecX,SpecY,groupnum);
    else
        set(handles.addpk,'enable','off');
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%
    % Artery---Vein----Separation
    if handles.enableAVS==1
        c=str2num(get(handles.SlideWindow,'string'));
        handles.GroupAVSVeloImg=zeros(OriFigX,OriFigY,groupnum);
    else
        set(handles.AVS,'enable','off');
    end
    %%%%%%%%%%

    %     set(handles.checkSV,'Enable','off');
    if handles.enableSK==1 && handles.enableSV==1
        handles.GroupSVeloImg=zeros(SpecX,SpecY,groupnum);
    else
        set(handles.SV,'Enable','off');
    end

    %     set(handles.checkTK,'Enable','off');
    if handles.enableTK==1
        handles.GroupTVeloImg=zeros(OriFigX,OriFigY,groupnum);
    else
        set(handles.TK,'Enable','off');
    end
    
        %     set(handles.checkSV,'Enable','off');
    if handles.enableTK==1 && handles.enableTV==1
        handles.GroupTtVeloImg=zeros(OriFigX,OriFigY,groupnum);
    else
        set(handles.TV,'Enable','off');
    end

    %     set(handles.checkSub,'Enable','Off');
    if handles.framenum>1 && handles.enableTSub==1
        handles.GroupTSubImg=zeros(OriFigX,OriFigY,groupnum);
    end


    StrTempFrame=fname((length(fname)-7):(length(fname)-4));
    StrTempGroup=fname(1:strfind(fname, intervChar));

    hwaitbar = waitbar(0,'Please wait...');
    tic;
    
    meanKvalue= zeros(SpecX,SpecY);
    meanaddvalue=zeros(SpecX,SpecY);
    meanNvalue= zeros(SpecX,SpecY);
    meansvalue= zeros(SpecX,SpecY);
    meandkvalue=zeros(SpecX,SpecY);
    SpeckleImg=zeros(OriFigX,OriFigY,handles.framenum);
    SpecklenImg=zeros(OriFigX,OriFigY,handles.framenum);
    SpeckledImg=zeros(OriFigX,OriFigY,handles.framenum);
    SpeckledkImg=zeros(OriFigX,OriFigY,handles.framenum);    
    for i=1:groupnum
        for j=1:framenum
            if loadmodelindex==1
                filename=[num2str(firstgroupindex+i-1),intervChar,...
                    sprintf(['%0',num2str(eval('filelen')),'d'],...
                    firstfileindex+j-1),fname(strfind(fname,'.'):length(fname))];
            elseif loadmodelindex==2
                filename=[num2str(firstgroupindex),intervChar,...
                    sprintf(['%0',num2str(eval('filelen')),'d'],...
                    (i-1)*framenum+firstfileindex+j-1),...
                    fname(strfind(fname,'.'):length(fname))];
            end
            fileread=[sprintf('%s',pname),filename];
            if strcmp(handles.load_filetype,'b16')
                SpeckleImg(:,:,j)=readb16(fileread);
                SpecklenImg(:,:,j)=readb16(fileread);
            else
                SpeckleImg(:,:,j)=imread(fileread);
                SpecklenImg(:,:,j)=imread(fileread);
            end
     
            if handles.enablenLASCA==1
                nn=3;
                % 标准化的散斑衬比计算方法，使用前3副散斑数据作为操作图像
                if j>=3
                    d=SpecklenImg(:,:,j-(nn-1))+SpecklenImg(:,:,j-(nn-2))+...
                    SpecklenImg(:,:,j-(nn-3)); 
                elseif j==1
                    d=3*SpecklenImg(:,:,1);
                elseif j==2
                    d=2*SpecklenImg(:,:,1)+SpecklenImg(:,:,2);
                end
                SpecklenImg(:,:,j)=d/3;
                an=spec2contr(double(SpecklenImg(:,:,j)),SlideWindow);
                meanNvalue=meanNvalue+an;
            end
            
            % 标准化的散斑衬比计算方法，对前（n）帧图片进行平均
            %if handles.enablenLASCA==1
                %d=0;
                %for k=1:j;
                    %d=SpecklenImg(:,:,j-k+1)+d;
                %end
                %SpecklenImg(:,:,j)=d/j;
                %an=spec2contr(double(SpecklenImg(:,:,j)),SlideWindow);
                %meanNvalue=meanNvalue+an;
            %end
            
            if handles.enableSK==1
                a=spec2contr(double(SpeckleImg(:,:,j)),SlideWindow);
                meanKvalue = meanKvalue+a;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % Pavel Zakharov and Andreas Volke,Alfred Buck,Bruno Weber,and Frank Scheffold,
            % Quantitative modeling of laser speckle imaging.Optics Letters 31(23):4265-4267 (2006)
            if handles.enableaddpt==1
                if j==1
                    P(:,:,j)=((0.5*(SpeckleImg(:,:,j)+SpeckleImg(:,:,j))).^2)./SpeckleImg(:,:,j)./SpeckleImg(:,:,j);
                elseif j==2
                    P(:,:,j)=((0.5*(SpeckleImg(:,:,j-1)+SpeckleImg(:,:,j))).^2)./SpeckleImg(:,:,j-1)./SpeckleImg(:,:,j);
                else
                    P(:,:,j)=((0.5*(SpeckleImg(:,:,j-2)+SpeckleImg(:,:,j))).^2)./SpeckleImg(:,:,j-2)./SpeckleImg(:,:,j);
                end
                P(:,:,j)=1./P(:,:,j);
                SpeckledImg(:,:,j)=SpeckleImg(:,:,j).*P(:,:,j);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % J.D.Briers.The statistics of fluctuating speckle patterns produced by a mixture of moving and stationary scatterers.
            % Optical and Quantum Electronics.10(4):364-366(1978)
            %空间衬比   
            if handles.enableaddpk==1
                n=handles.framenum;
                dk=sqrt(sum(SpeckleImg.^2,3)/(n-1)-n*mean(SpeckleImg,3).^2/(n-1))./mean(SpeckleImg,3);
                % 使用空间衬比计算
                %dk(:,:,j)=spec2contr(double(SpeckleImg(:,:,j)),SlideWindow);
                %ddk=mean(dk,3);
                Q=1-sqrt(abs(1-dk));
                %h=str2num(get(handles.SlideWindow,'string'));
                %SpeckledkImg((h+1)/2:OriFigX-(h-1)/2,(h+1)/2:OriFigY-(h-1)/2,j)=SpeckleImg((h+1)/2:OriFigX-(h-1)/2,(h+1)/2:OriFigY-(h-1)/2,j).*Q;
                % 使用时间衬比计算
                SpeckledkImg(:,:,j)=SpeckleImg(:,:,j).*Q;
                % 使用前几幅图片的均值作为操作图像
                %nnn=3;
                %if j>=3
                %    dd=SpeckledkImg(:,:,j-(nnn-1))+SpeckledkImg(:,:,j-(nnn-2))+SpeckledkImg(:,:,j-(nnn-3)); 
                %elseif j==1
                %    dd=3*SpecklenImg(:,:,1);
                %elseif j==2
                %    dd=2*SpecklenImg(:,:,1)+SpecklenImg(:,:,2);
                %end
                %SpeckledkImg(:,:,j)=dd/3;                                
                ck=spec2contr(double(SpeckledkImg(:,:,j)),SlideWindow);
                meandkvalue=meandkvalue+ck;
            end
                           
            if handles.enablesLASCA==1
                b=spec2contr(double(SpeckleImg(:,:,j)),SlideWindow);
                meansvalue=meansvalue+b;
            end
            
            
            % 时间+空间联合衬比计算方法：Slidewindow*Slidewindow*3(标准模式，但是数据量巨大，out of memory)
            % 空间衬比计算方法与时间衬比计算方法的误差达到13%，
            % 而时空联合衬比计算方法的误差只有5%
            % if handles.enablestLASCA==1
            % c=str2num(get(handles.SlideWindow,'string'));
            % d=3;
            % p=zeros(c,c,d);
            % C=zeros(c*c*d,1);
            % m=OriFigX;
            % n=OriFigY;
            %pq=zeros(OriFigX,OriFigY,d);
            %if j<=d-1
            %    pq=SpeckleImg(:,:,j:j+d-1);
            %else
            %    pq=SpeckleImg(:,:,j-d+1:j);
            %end
            %for h=((c+1)/2):(m-(c-1)/2)
            %    for k=((c+1)/2):(n-(c-1)/2)
            %        p=pq(h-(c-1)/2:h+(c-1)/2,k-(c-1)/2:k+(c-1)/2,:);
            %        for q=1:d
            %            A=p(:,:,d);
            %            B=reshape(A,c*c,1);
            %            C((q-1)*c*c+1:q*c*c)=B;
            %        end
            %        D=mean(C);
            %        E=std(C);
            %        handles.GroupstVeloImg(h-(c-1)/2,k-(c-1)/2,i)=E/D;
            %    end
            %end
            %handles.GroupstVeloImg=handles.GroupstVeloImg(1:m-(c-1),1:n-(c-1),:)
            %end          
        end
                                                  
        if handles.enablestLASCA==1
            c=str2num(get(handles.SlideWindow,'string'));
            p=zeros(c,c,handles.framenum);
            C=zeros(c*c*handles.framenum,1);
            m=OriFigX;
            n=OriFigY;
            for h=((c+1)/2):(m-(c-1)/2)
                for k=((c+1)/2):(n-(c-1)/2)
                    p=SpeckleImg(h-(c-1)/2:h+(c-1)/2,k-(c-1)/2:k+(c-1)/2,:);
                    for q=1:handles.framenum
                        A=p(:,:,q);
                        B=reshape(A,c*c,1);
                        C((q-1)*c*c+1:q*c*c)=B;
                    end
                    D=mean(C);
                    E=std(C);
                    handles.GroupstVeloImg(h-(c-1)/2,k-(c-1)/2,i)=E/D;
                end
            end
            handles.GroupstVeloImg=handles.GroupstVeloImg(1:m-(c-1),1:n-(c-1),:);
        end      
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        if handles.enableAVS==1
            PI=3.1416;
            c=str2num(get(handles.SlideWindow,'string'));
            n=handles.framenum;
            mm=OriFigX;
            nn=OriFigY;
            av=zeros(1,handles.framenum);
            B=zeros(OriFigX,OriFigY);
            D=zeros(OriFigX,OriFigY);
            E=zeros(OriFigX,OriFigY);
            G=zeros(10,15);
            H=zeros(5,5);
            a=0;
            b=0;
            aa=0;
            bb=0;
            
            %%%  step No.1
            for hh=1:OriFigX;
                for kk=1:OriFigY
                    av=SpeckleImg(hh,kk,:);
                    B(hh,kk)=min(av);%实际的激光散斑最小光强
                    %B(hh,kk)=mean(av)-std(av)*sqrt(PI/(4-PI));%理论上的激光散斑最小光强
                end
            end
            
            %%%  step No.2
            A=sqrt(sum(SpeckleImg.^2,3)/(n-1)-n*mean(SpeckleImg,3).^2/(n-1))./mean(SpeckleImg,3);
            M=A;
            D=mean(SpeckleImg,3); 
            T=0.1;% 设定阈值，将血管与背景分开
            for h=1:OriFigX
                for k=1:OriFigY
                    if A(h,k)>=T %大于阈值，非血管区域，用“1”表示；
                        A(h,k)=1;
                    else
                        A(h,k)=0;%小于阈值，血管区域,用“0”表示；
                    end
                end
            end
            
            D=D.*A;
            for h=1:OriFigX-50
                for k=1:OriFigY
                    if k<=500
                        G=D(h:h+9,k:k+14);
                    else
                        G=D(h:h+9,k-14:k);
                    end
                    if A(h,k)==0
                        for p=1:10
                            for q=1:15
                                if G(p,q)~=0
                                    a=a+1;
                                    b=b+G(p,q);
                                end
                            end
                        end
                        E(h,k)=b/a;
                    else
                        E(h,k)=D(h,k);% computing the background image.
                    end
                    a=0;
                    b=0;
                end
            end
            
            for h=OriFigX-49:OriFigX
                for k=1:OriFigY
                    if k<=500
                        G=D(h-9:h,k:k+14);
                    else
                        G=D(h-9:h,k-14:k);
                    end                           
                    if A(h,k)==0
                        for p=1:10
                            for q=1:15
                                if G(p,q)~=0
                                    a=a+1;
                                    b=b+G(p,q);
                                end
                            end
                        end
                        E(h,k)=b/a;
                        %E(h-(c-1)/2,k-(c-1)/2)=mean(mean(D(h-(c-1)/2:h+(c-1)/2,k-(c-1)/2:k+(c-1)/2)));% computing the background image.
                    else
                        E(h,k)=D(h,k);% computing the background image.
                    end
                    a=0;
                    b=0;
                end
            end
            
            F=E;
            %d=zeros(5,5);
            %for i=3:518
            %    for j=3:694
            %        d=E(i-2:i+2,j-2:j+2);
            %        F(i,j)=mean(mean(d));
            %    end
            %end 
            MR=B./F;%computing the relative temporal minimum reflectance image.
            RMR=B./F.*A;

            
            % step No.3 :Artery_Vein separation.computing the RTMR
            for h=1:OriFigX-50
                for k=1:OriFigY
                    if k<=500
                        H=RMR(h:h+4,k:k+4);
                    else
                        H=RMR(h:h+4,k-4:k);
                    end
                    if A(h,k)==0
                        for p=1:5
                            for q=1:5
                                if H(p,q)~=0
                                    aa=aa+1;
                                    bb=bb+H(p,q);
                                end
                            end
                        end
                        RTMR(h,k)=bb/aa;
                    else
                        RTMR(h,k)=MR(h,k);% computing the background image.
                    end
                    aa=0;
                    bb=0;
                end
            end
            
            for h=OriFigX-49:OriFigX
                for k=1:OriFigY
                    if k<=500
                        H=RMR(h-4:h,k:k+4);
                    else
                        H=RMR(h-4:h,k-4:k);
                    end                           
                    if A(h,k)==0
                        for p=1:5
                            for q=1:5
                                if H(p,q)~=0
                                    aa=aa+1;
                                    bb=bb+H(p,q);
                                end
                            end
                        end
                        RTMR(h,k)=bb/aa;
                    else
                        RTMR(h,k)=MR(h,k);% computing the background image.
                    end
                    aa=0;
                    bb=0;
                end
            end
            
            %step No.4: Artery_Vein separation.
            for h=1:OriFigX
                for k=1:OriFigY
                    if A(h,k)==0
                    if RTMR(h,k)<=MR(h,k)
                        R(h,k)=1;%动脉点
                    else
                        R(h,k)=0;%静脉点
                    end
                    if R(h,k)==1
                        M(h,k)=0.8;
                    else
                        M(h,k)=0;
                    end
                    else
                        M(h,k)=M(h,k);
                    end
                end
            end
            N=M;
                                   
            % 判断这一点的划分是否正确，如果这点的3*3区域里有n(n>=3)以上为1
            % 则这点的划分正确，若(n<3)说明这点的划分错误
            %cc=0;
            %d=zeros(3,3);
            %for i=2:OriFigX-1
            %    for j=2:OriFigY-1
            %        if M(i,j)==0.6
            %            d=M(i-1:i+1,j-1:j+1);
            %        end
            %        for p=1:3
            %            for q=1:3
            %                if d(p,q)==0.6
            %                    cc=cc+1;
            %                end
            %            end
            %        end
            %        if cc<=3
            %            N(i,j)=mean(mean(d));
            %        end
            %        cc=0;
            %    end
            %end
            handles.GroupAVSVeloImg(:,:,i)=N;
        end
            %%%%%%%%%%%%%%%%%%%%%%%% 
            
        %clear tmp a
        handles.GroupSpecImg(:,:,i)=mean(SpeckleImg,3);

        if handles.enableSK==1
            meanKvalue=meanKvalue/handles.framenum;
            handles.GroupKImg(:,:,i)=meanKvalue;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.enableaddpt==1
            n=handles.framenum;
            handles.GroupaddVeloImg(:,:,i)=sqrt(sum(SpeckledImg.^2,3)/(n-1)-n*mean(SpeckledImg,3).^2/(n-1))./mean(SpeckledImg,3);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if handles.enableaddpk==1
        % 归一化处理
        % meandkvalue=(meandkvalue-min(min(meandkvalue)))./(max(max(meandkvalue))-min(min(meandkvalue)));
            meandkvalue=meandkvalue/handles.framenum;
            handles.GroupaddpkImg(:,:,i)=meandkvalue;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        if handles.enablenLASCA==1
            meanNvalue=meanNvalue/handles.framenum;
            handles.GroupnImg(:,:,i)=meanNvalue;
        end
        
        if handles.enablesLASCA==1
            meansvalue=meansvalue/handles.framenum;
            handles.GroupsImg(:,:,i)=meansvalue;
        end

        if handles.enableSK==1 && handles.enableSV==1
           handles.GroupSVeloImg(:,:,i)=velocity(SpecX*SpecY,meanKvalue);
        end

        if handles.enableTK==1
            n=handles.framenum;
            %SpeckleImg2=SpeckleImg.^2;
            %MeanSpeckleImg=mean(SpeckleImg,3);
            %MeanSpeckleImg2=mean(SpeckleImg2,3);
            %handles.GroupTVeloImg(:,:,i)=sqrt(n*MeanSpeckleImg2/(n-1)-n*MeanSpeckleImg.^2/(n-1))./MeanSpeckleImg;
            TtV=sqrt(sum(SpeckleImg.^2,3)/(n-1)-n*mean(SpeckleImg,3).^2/(n-1))./mean(SpeckleImg,3);
            handles.GroupTVeloImg(:,:,i)=sqrt(sum(SpeckleImg.^2,3)/(n-1)-n*mean(SpeckleImg,3).^2/(n-1))./mean(SpeckleImg,3);
        end
        
        if handles.enableTK==1 && handles.enableTV==1
           handles.GroupTtVeloImg(:,:,i)=velocity(OriFigX*OriFigY,TtV);
        end
                

        if handles.enabletLASCA==1;
            n=handles.framenum;
            handles.GrouptVeloImg(:,:,i)=sqrt(sum(SpeckleImg.^2,3)/(n-1)-n*mean(SpeckleImg,3).^2/(n-1))./mean(SpeckleImg,3);
        end
        
        %keyboard;

        if handles.framenum>1 && handles.enableTSub==1
            numerator=0;
            denominator=0;
            for subi=2:framenum
                %handles.GroupTSubImg(:,:,i)=handles.GroupTSubImg(:,:,i)+abs(SpeckleImg(:,:,subi)-SpeckleImg(:,:,subi-1))./abs(SpeckleImg(:,:,subi)+SpeckleImg(:,:,subi-1));
                numerator=numerator+abs(SpeckleImg(:,:,subi)-SpeckleImg(:,:,subi-1));
                denominator=denominator+abs(SpeckleImg(:,:,subi)+SpeckleImg(:,:,subi-1));
            end
            %handles.GroupTSubImg=handles.GroupTSubImg./(framenum-1);
            handles.GroupTSubImg(:,:,i)=numerator./denominator;
        end

        %handles.GroupTVeloImg(:,:,i)=sqrt(MeanSpeckleImg.^2./(MeanSpeckleImg2-MeanSpeckleImg.^2+eps));
        
        t=toc;
        %keyboard
        meanKvalue=meanKvalue.*0;
        SpeckleImg=SpeckleImg.*0;
        ts=(groupnum-i)*t/i;
        str = sprintf('%d %s %d %s',floor(ts/60), 'minutes and ', floor(mod(ts,60)),' seconds left...');
        waitbar(i/groupnum, hwaitbar, str);
    end
        
    close(hwaitbar);
    handles.cur_data=handles.GroupSpecImg(:,:,1);
    
    imagesc(handles.cur_data);
    
    %     keyboard
    
    axis off;
    axis equal;
    colormap(gray)
    lim=min(min(handles.cur_data));
    uim=max(max(handles.cur_data));
    caxis([lim,uim]);
    handles.pname=pname;
    handles.fname=fname;
    handles.Cimage=getimage;
    handles.SpecRectX=SpecX;
    handles.SpecRectY=SpecY;
    handles.OrigRectX=OriFigX;
    handles.OrigRectY=OriFigY;

    [m,n]=size(handles.cur_data);
    frame_width=n;
    frame_height=m;
    handles.frame_width=frame_width;
    handles.frame_height=frame_height;
   
    set(handles.figure1,'name',[sprintf('Speckle statistic analysis--'),pname]);
    set(handles.figure1,'WindowButtonMotionFcn','SpecUI(''NormalMotionFcn'',gcbo,[],guidata(gcbo))');
   
    % set(handles.figure1,'WindowButtonMotionFcn','SpecUI(''NormalMotionFcn'')');
    % set(gcf,'WindowButtonMotionFcn','NormalMotionFcn');
    
    guidata(gcbo,handles);
       
    set(hObject,'Enable','Off');
else
    errordlg('Please set the frame and group number','Error!')
end

function Speckle_Callback(hObject, eventdata, handles)
handles.cur_data = handles.GroupSpecImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['speckle image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;
handles.VImg=0;
guidata(gcbo,handles);
axis off;
axis equal;

% --- Executes on button press in checkSK.
function checkSK_Callback(hObject, eventdata, handles)
handles.enableSK=~handles.enableSK;

if handles.enableSK==1
    set(handles.SK,'enable','on');
else
    set(handles.SK,'enable','off');
    set(handles.checkSV,'enable','on','value',0);
    set(handles.SV,'enable','off');
end
guidata(gcbo,handles);

function SK_Callback(hObject, eventdata, handles)
handles.cur_data = handles.GroupKImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['spatial speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;
handles.VImg=1;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;
% keyboard

% --- Executes on button press in checkSV.
function checkSV_Callback(hObject, eventdata, handles)
handles.enableSV=~handles.enableSV;

if handles.enableSV==1
    set(handles.SK,'enable','on');
    handles.enableSK=1
    set(handles.checkSK,'enable','on','value',1);
    set(handles.SV,'enable','on');
else
    set(handles.SV,'enable','off');
end
guidata(gcbo,handles);

function SV_Callback(hObject, eventdata, handles)
handles.cur_data = handles.GroupSVeloImg(:,:,handles.currentgroup);
imagesc(handles.cur_data);colormap(gray),title(['blood velocity image by contrast: group ',...
    sprintf('%d',handles.currentgroup)]);
handles.VImg=2;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;

% --- Executes on button press in checkTK.
function checkTK_Callback(hObject, eventdata, handles)
handles.enableTK=~handles.enableTK;

if handles.enableTK==1
    set(handles.SK,'enable','on');
else
    set(handles.TK,'enable','off');
    set(handles.checkTV,'enable','on','value',0);
    set(handles.TV,'enable','off');
end
guidata(gcbo,handles);


function TK_Callback(hObject, eventdata, handles)
handles.cur_data=handles.GroupTVeloImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['temporal speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]);colormap(gray);
handles.Cimage=getimage;
handles.VImg=3;
guidata(gcbo,handles);
axis off;
axis equal;
% keyboard

% --- Executes on button press in checkTV.
function checkTV_Callback(hObject, eventdata, handles)
handles.enableTV=~handles.enableTV;

if handles.enableTV==1
    set(handles.TK,'enable','on');
    handles.enableTK=1
    set(handles.checkTK,'enable','on','value',1);
    set(handles.TV,'enable','on');
else
    set(handles.TV,'enable','off');
end
guidata(gcbo,handles);

function TV_Callback(hObject, eventdata, handles)
handles.cur_data = handles.GroupTtVeloImg(:,:,handles.currentgroup);
imagesc(handles.cur_data);colormap(gray),title(['blood velocity image by contrast: group ',...
    sprintf('%d',handles.currentgroup)]);

%%%%%%%
caxis([0 400]); %调节TV显示的colorbar,但并不会对TV的真实值造成影响！
%%%%%%%

handles.VImg=4;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;

% --- Executes on button press in checkSub.
function checkSub_Callback(hObject, eventdata, handles)
handles.enableTSub=~handles.enableTSub;

if handles.enableTSub==1
    set(handles.Sub,'enable','on');
else
    set(handles.Sub,'enable','off');
end
guidata(gcbo,handles);

function Sub_Callback(hObject, eventdata, handles)
handles.cur_data=handles.GroupTSubImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['Statistic imformation by subtraction ',...
    sprintf('%d',handles.currentgroup)]);colormap(gray);
handles.Cimage=getimage;
guidata(gcbo,handles);
axis off;
axis equal;


function SetROI_Callback(hObject, eventdata, handles)
handles.rect = ceil(getrect(gcf));
if handles.VImg==0
    [handles.OrigRectX handles.OrigRectY]=size(imcrop(getimage,handles.rect));
    imagesc(imcrop(getimage,handles.rect));colormap(gray);
elseif handles.VImg==1
    [handles.SpecRectX handles.SpecRectY]=size(imcrop(getimage,handles.rect));
    imagesc(imcrop(getimage,handles.rect));colormap(gray);
elseif handles.VImg==2
    [handles.SpecRectX handles.SpecRectY]=size(imcrop(getimage,handles.rect));
    imagesc(imcrop(getimage,handles.rect));colormap(gray);
elseif handles.VImg==3
    [handles.OrigRectX handles.OrigRectY]=size(imcrop(getimage,handles.rect));
    imagesc(imcrop(getimage,handles.rect));colormap(gray);
elseif handles.VImg==4
    [handles.OrigRectX handles.OrigRectY]=size(imcrop(getimage,handles.rect));
    imagesc(imcrop(getimage,handles.rect));colormap(gray);
end
handles.Cimage=getimage;
guidata(gcbo,handles);

%skip
function Frames_Callback(hObject, eventdata, handles)
function Frames_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
function Groups_Callback(hObject, eventdata, handles)
function Groups_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Save_Callback(hObject, eventdata, handles)
[savename, savepath]= uiputfile({'*.mat';'*.png';'*.tif';'*.jpg';'*.gif'},'Save current figure!');
if (length(savename)==0)|(length(savepath)==0)
    return;
else
    filedetect = findstr(savename,'.');
    fmt = savename((filedetect+1):length(savename));
    A = getimage(gcf);
    mintmp=min(min(A));
    maxtmp=max(max(A));
    tmp = (A-mintmp)*1.0/(maxtmp-mintmp); % 归一化处理，量化基线！
    m = uint16(round(4096*tmp));
    n = uint8(round(255*tmp));
    file = [savepath,savename];
    switch (fmt)
        case 'mat'
            save(savename,'A');
        case 'png'
            imwrite(m, file, 'BitDepth', 16);
        case 'tif'
            imwrite(m, file, 'tif');
        case 'jpg'
            cmap=colormap;
            imwrite(n,file,'jpg');
        case 'gif'
            imwrite(n, file, 'gif');
    end
end

% keyboard
% save savename A;
% function Vessel_Callback(hObject, eventdata, handles)
% if handles.setroi==1
%     [cx,cy,c]=improfile;
%     cx=ceil(cx);
%     cy=ceil(cy);
%     DynaV=zeros(handles.groupnum,length(cx));
%     if handles.algorithm==1
%         for i=1:handles.groupnum
%             for k=1:length(cx)
%                 DynaV(i,k)=velo(cx(k), cy(k));
%             end
%         end
%         imagesc(DynaV);
%     elseif handles.algorithm==0
%          for i=1:handles.groupnum
%             for k=1:length(cx)
%                 DynaV(i,k)=S(cx(k), cy(k));
%             end
%         end
%         imagesc(DynaV);
%     else
%         errordlg('Please show the velocity image','Serious Problem')
%     end
% else
%     errordlg('You did not set ROI','Serious Problem')
% end

function Grid_Callback(hObject, eventdata, handles)
grid
grid minor

function OrgImg_CreateFcn(hObject, eventdata, handles)
function pushbutton11_Callback(hObject, eventdata, handles)
function pushbutton12_Callback(hObject, eventdata, handles)

% --- Executes on slider movement.
function BlackThreshold_Callback(hObject, eventdata, handles)
NewVal=get(hObject,'value');
set(handles.BlackThreshold,'string',round(NewVal));
% keyboard;
%
maxval=max(max(handles.Cimage));
minval=min(min(handles.Cimage));
A=ceil(((handles.Cimage-minval)*1.0/(maxval-minval))*4096);
imagesc(A)
%
cmax=str2double(get(handles.WhiteThreshold,'string'));
cmin=str2double(get(handles.BlackThreshold,'string'));
if cmax>cmin
    caxis([cmin cmax]);
end
guidata(gcbo,handles);
axis off;
axis equal;

function BlackThreshold_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function WhiteThreshold_Callback(hObject, eventdata, handles)
NewVal=get(hObject,'value');
set(handles.WhiteThreshold,'string',round(NewVal));
%
maxval=max(max(handles.Cimage));
minval=min(min(handles.Cimage));
A=ceil((handles.Cimage-minval)*4096/(maxval-minval));
imagesc(A)
%
cmax=str2double(get(handles.WhiteThreshold,'string'));
cmin=str2double(get(handles.BlackThreshold,'string'));
if cmax>cmin
    caxis([cmin cmax]);
end
guidata(gcbo,handles);
axis off;
axis equal;


function WhiteThreshold_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Gray2Color_Callback(hObject, eventdata, handles)
contents=get(hObject,'value');
switch contents
    case 1.0
        imagesc(handles.Cimage);
        colormap(gray);
    otherwise
        imagesc(handles.Cimage);
        colormap(jet(64));
end
axis off;
axis equal;
guidata(gcbo,handles);

function Gray2Color_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BW_Callback(hObject, eventdata, handles)
imagesc(handles.Cimage);
if (get(hObject,'Value') == get(hObject,'Max'))
    I=handles.Cimage;
    maxval=max(max(handles.Cimage));
    minval=min(min(handles.Cimage));
    A=ceil((handles.Cimage-minval)*4096/(maxval-minval));
    level=300;
    bw = A>=level;
    imagesc(bw)
    set(handles.BWThreshold,'Enable','On');
    set(handles.BWThreshold,'Value',level/4096);
else
    imagesc(handles.Cimage);
    set(handles.BWThreshold,'Enable','Off');
end
guidata(gcbo,handles);
axis off;
axis equal;

function BWThreshold_Callback(hObject, eventdata, handles)
NewVal=get(hObject,'value');
set(handles.BWThreshold,'string',NewVal);
imagesc(handles.Cimage);
I=handles.Cimage;
maxval=max(max(handles.Cimage));
minval=min(min(handles.Cimage));
A=ceil((handles.Cimage-minval)*4096/(maxval-minval));
%     bw = im2bin(A,4096*NewVal);
level=4096*NewVal
bw = A>=level;
imagesc(bw)
axis off;
axis equal;

function BWThreshold_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function LoadEnable_Callback(hObject, eventdata, handles)
set(handles.Load,'Enable','On');


% --- Executes on button press in speckle.
function speckle_Callback(hObject, eventdata, handles)
% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function MultiROI_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
    handles.bmroi=1;
else
    handles.bmroi=0;
end
guidata(gcbo,handles);

function LoadModel_Callback(hObject, eventdata, handles)
set(handles.LoadEnable,'Enable','on');
% contents=get(hObject,'value');
% switch contents
%     case 1.0
%         imagesc(handles.Cimage);
%         colormap(gray);
%     otherwise
%         imagesc(handles.Cimage);
%         colormap(jet(64));
% end
% guidata(gcbo,handles);

function LoadModel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function First_Callback(hObject, eventdata, handles)
handles.currentgroup=handles.FirstGroup;
guidata(gcbo,handles);
if handles.VImg==0
    SpecUI('Speckle_Callback',gcbo,[],guidata(gcbo))
    %     Speckle_Callback(gcbo,[],guidata(gcbo))
elseif handles.VImg==1
    SpecUI('SK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==2
    SpecUI('SV_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==3
    SpecUI('TK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==4
    SpecUI('TV_Callback',gcbo,[],guidata(gcbo));
end


function Rewind_Callback(hObject, eventdata, handles)
handles.currentgroup=handles.currentgroup-1;
if handles.currentgroup<1;
    handles.currentgroup=1;
end
guidata(gcbo,handles);
if handles.VImg==0
    % SpecUI('Speckle_Callback',gcbo,[],guidata(gcbo))
    Speckle_Callback(gcbo,[],guidata(gcbo))
elseif handles.VImg==1
    SpecUI('SK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==2
    SpecUI('SV_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==3
    SpecUI('TK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==4
    SpecUI('TV_Callback',gcbo,[],guidata(gcbo));
end

function Forward_Callback(hObject, eventdata, handles)
handles.currentgroup=handles.currentgroup+1;
if handles.currentgroup>handles.LastGroup;
    handles.currentgroup=handles.LastGroup;
end
guidata(gcbo,handles);
if handles.VImg==0
    % SpecUI('Speckle_Callback',gcbo,[],guidata(gcbo))
    Speckle_Callback(gcbo,[],guidata(gcbo))
elseif handles.VImg==1
    SpecUI('SK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==2
    SpecUI('SV_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==3
    SpecUI('TK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==4
    SpecUI('TV_Callback',gcbo,[],guidata(gcbo));

end

function Last_Callback(hObject, eventdata, handles)
handles.currentgroup=handles.LastGroup;
guidata(gcbo,handles);
if handles.VImg==0
    % SpecUI('Speckle_Callback',gcbo,[],guidata(gcbo))
    Speckle_Callback(gcbo,[],guidata(gcbo))
elseif handles.VImg==1
    SpecUI('SK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==2
    SpecUI('SV_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==3
    SpecUI('TK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==4
    SpecUI('TV_Callback',gcbo,[],guidata(gcbo));
end

function Current_Callback(hObject, eventdata, handles)
NewVal=get(hObject,'value');
% set(handles.Current,'string',round(NewVal));
handles.currentgroup=round(NewVal);
if handles.currentgroup>handles.LastGroup;
    handles.currentgroup=handles.LastGroup;
elseif handles.currentgroup<1;
    handles.currentgroup=1;
end
set(hObject,'Value',handles.currentgroup);
set(handles.EditCurrent,'String',num2str(handles.currentgroup));
guidata(gcbo,handles);
if handles.VImg==0
    % SpecUI('Speckle_Callback',gcbo,[],guidata(gcbo))
    Speckle_Callback(gcbo,[],guidata(gcbo))
elseif handles.VImg==1
    SpecUI('SK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==2
    SpecUI('SV_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==3
    SpecUI('TK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==4
    SpecUI('TV_Callback',gcbo,[],guidata(gcbo));
end


function Current_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function EditCurrent_Callback(hObject, eventdata, handles)
NewVal=get(hObject,'string');
if length(NewVal) == 0
    handles.currentgroup=1;
else
    handles.currentgroup=str2num(NewVal);
    if handles.currentgroup <1
        handles.currentgroup=1;
    elseif handles.currentgroup>handles.LastGroup
        handles.currentgroup=handles.LastGroup;
    end
end
guidata(gcbo,handles);
if handles.VImg==0
    % SpecUI('Speckle_Callback',gcbo,[],guidata(gcbo))
    Speckle_Callback(gcbo,[],guidata(gcbo))
elseif handles.VImg==1
    SpecUI('SK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==2
    SpecUI('SV_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==3
    SpecUI('TK_Callback',gcbo,[],guidata(gcbo));
elseif handles.VImg==4
    SpecUI('TV_Callback',gcbo,[],guidata(gcbo));
end

function EditCurrent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Plot_Callback(hObject, eventdata, handles)
handles.originhandle=gcf;
handles.plothandle=figure;
figure(handles.originhandle)
i=1;
% if handles.bmroi == 1
%     while 1
%         w = waitforbuttonpress;
%         if w == 1
%             figure(handles.originhandle);
%             break;
%         else
%             rect(i,:)=ceil(getrect(gcf));
%             rectangle('Position',rect(i,:),...
%                 'LineWidth',1,'LineStyle',':','EdgeColor','r');
%             text(rect(i,1),rect(i,2)-10,num2str(i),'color','r');
%             i=i+1;
%         end
%     end
%     rectnum=i-1;
if handles.bmroi == 1
    while 1
        rect(i,:)=ceil(getrect(gcf));
        rectangle('Position',rect(i,:),...
            'LineWidth',1,'LineStyle',':','EdgeColor','r');
        text(rect(i,1),rect(i,2)-10,num2str(i),'color','r');
        i=i+1;
        w = waitforbuttonpress;
        if w == 1
            figure(handles.originhandle);
            break;
        end
    end
    rectnum=i-1;
else
    rect(i,:)=ceil(getrect(gcf));
    rectangle('Position',rect(i,:),...
        'LineWidth',1,'LineStyle',':','EdgeColor','r');
    text(rect(i,1),rect(i,2)-10,num2str(i),'color','r');
    rectnum=i;
end

DynaV=zeros(handles.groupnum,rectnum);
for m=1:rectnum
    if handles.VImg==0
        for i=1:handles.groupnum
            DynaV(i,m)=mean2(imcrop(handles.GroupSpecImg(:,:,i),rect(m,:)));
        end
    elseif handles.VImg==1
        for i=1:handles.groupnum
            DynaV(i,m)=mean2(imcrop(handles.GroupKImg(:,:,i),rect(m,:)));
        end
        a=1
    elseif handles.VImg==2
        for i=1:handles.groupnum
            DynaV(i,m)=mean2(imcrop(handles.GroupSVeloImg(:,:,i),rect(m,:)));
        end
    elseif  handles.VImg==3
        for i=1:handles.groupnum
            DynaV(i,m)=mean2(imcrop(handles.GroupTVeloImg(:,:,i),rect(m,:)));
        end
        a=2
    elseif  handles.VImg==4
        for i=1:handles.groupnum
            DynaV(i,m)=mean2(imcrop(handles.GroupTtVeloImg(:,:,i),rect(m,:)));
        end
        guidata(gcbo,handles);
    else
        errordlg('Please show the contrast or velocity image','Serious Problem')
    end
end

figure(handles.plothandle);
plot(DynaV);
xlabel('group index');
ylabel('relative value');
if handles.VImg==0
    title('Raw Speckle Dynamics');
elseif  handles.VImg==1
    title('Spatial Statistics Contrast Dynamics');
elseif  handles.VImg==2
    title('Spatial Statistics Velocity Dynamics');
elseif  handles.VImg==3
    title('Temporal Statistics Contrast Dynamics');
elseif  handles.VImg==4
    title('Temporal Statistics Velocity Dynamics');
end
legend('show');
figure(handles.originhandle);


% if handles.VImg==0|handles.VImg==1|handles.VImg==2|handles.VImg==3
%     button = questdlg('Save the lines?','Save lines','Yes','No','Yes');
%     if length(button) == length('Yes')
%         [filename, pathname] = uiputfile('*.mat', 'Save the lines');
%         if isequal(filename,0) | isequal(pathname,0)
%             disp('User selected Cancel')
%         else
%             index=findstr(filename,'.');
%             for i=1:rectnum
%                 fileid=[filename(1:index-1),num2str(i),filename(index:length(filename))];
%                 keyboard;
%                 m=DynaV(i,:)
%                 eval(['save ' fileid ' m   -ascii']);
%             end
%         end
%     else return;
%     end
% end



% --- Executes on button press in ST.
function ST_Callback(hObject, eventdata, handles)
handles.originhandle=gcf;
handles.plothandle=figure;
figure(handles.originhandle)
i=1;

if handles.bmroi == 1
    while 1
        rect(i,:)=ceil(getrect(gcf));
        rectangle('Position',rect(i,:),...
            'LineWidth',1,'LineStyle',':','EdgeColor','r');
        text(rect(i,1),rect(i,2)-10,num2str(i),'color','r');
        i=i+1;
        w = waitforbuttonpress;
        if w == 1
            figure(handles.originhandle);
            break;
        end
    end
    rectnum=i-1;
else
    rect(i,:)=ceil(getrect(gcf));
    rectangle('Position',rect(i,:),...
        'LineWidth',1,'LineStyle',':','EdgeColor','r');
    text(rect(i,1),rect(i,2)-10,num2str(i),'color','r');
    rectnum=i;
end

DynaV=zeros(handles.groupnum,rectnum);
if handles.VImg==3
    for m=1:rectnum
        Temprect = zeros(size(imcrop(handles.GroupTVeloImg(:,:,i),rect(m,:))));
        %         keyboard
        for i=1:handles.groupnum
            Temprect=Temprect+imcrop(handles.GroupTVeloImg(:,:,i),rect(m,:));
        end

        index=find(Temprect>mean2(Temprect));
        %         mean2(Temprect(index));

        for i=1:handles.groupnum
            temp=imcrop(handles.GroupTVeloImg(:,:,i),rect(m,:));
            DynaV(i,m)=mean2(temp(index));
        end
    end
end

figure(handles.plothandle);
plot(DynaV);
xlabel('Relative optical intensity');
ylabel('Temporal speckle contrast');
figure(handles.originhandle);


% --- Executes on button press in PC.
function PC_Callback(hObject, eventdata, handles)
handles.originhandle=gcf;
handles.plothandle=figure;
figure(handles.originhandle)
tmp=ceil(ginput(1));
x=tmp(1);
y=tmp(2);
DynaV=zeros(1,handles.groupnum);
DynaVC=zeros(1,handles.groupnum);


if handles.VImg==0
    for i=1:handles.groupnum
        DynaV(i)=handles.GroupSpecImg(x,y,i);
    end
elseif handles.VImg==1
    for i=1:handles.groupnum
        DynaV(i)=handles.GroupKImg(x,y,i);
    end
elseif handles.VImg==2
    for i=1:handles.groupnum
        DynaV(i)=handles.GroupSVeloImg(x,y,i);
    end
elseif  handles.VImg==3
    for i=1:handles.groupnum
        DynaV(i)= handles.GroupTVeloImg(x,y,i) ;
    end
elseif  handles.VImg==4
    for i=1:handles.groupnum
        handles.GroupTtVeloImg(x,y,i);
    end
else
    errordlg('Please show the right image','Serious Problem')
end
figure(handles.plothandle);
plot(1:handles.groupnum,DynaV,'*k:');
xlabel('group index');
ylabel('relative value');
if handles.VImg==0
    title('Raw Speckle Dynamics');
elseif  handles.VImg==1
    title('Spatial speckle contrast dynamics');
elseif  handles.VImg==2
    title('Spatial speckle velocity dynamics');
elseif  handles.VImg==3
    title('Temporal speckle contrast dynamics');
elseif  handles.VImg==4
    title('Temporal speckle velocity dynamics');
end
legend('show');
figure(handles.originhandle);


% function NormalMotionFcn(arg)
% % Set the cursor to a Cross-Hair when above the Original image, and back
% % to an arrow when not.
% % This is the normal motion function for the window when we are not in
% % a MyGetline selection state.
% if ~nargin
%     handles=guidata(gcbo);
% %     pos=get(handles.OrgImg,'Position');
%     pt=get(handles.figure1,'CurrentPoint');
% pt=get(gcf,'CurrentPoint');
%     x=pt(1,1);
%     y=pt(1,2);
% 
%     if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4) )
%         cp=get(handles.OrgImg,'CurrentPoint');
% 
%         if round(cp(1,2))> 0 & round(cp(1,2))<handles.frame_height & round(cp(1,1))>0 & round(cp(1,1))<handles.frame_width
%             set(handles.figure1,'pointer','cross');
%             set(handles.mx,'string',cp(1,1));
%             set(handles.my,'string',cp(1,2));
%             set(handles.mvalue,'string',handles.cur_data(round(cp(1,2)),round(cp(1,1))));
%         else
%             set(handles.figure1,'pointer','arrow');
%         end
% 
%     else
%         set(handles.figure1,'pointer','arrow');
%     end
% elseif strcmp(arg,'getline')
%     handles=guidata(gcbo);
%     pos=get(handles.OrgImg,'Position');
%     pt=get(handles.figure1,'CurrentPoint');
%     x=pt(1,1);
%     y=pt(1,2);
%     if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4))
%         set(handles.figure1,'pointer','crosshair');
%         cp=get(handles.OrgImg,'CurrentPoint');
%         set(handles.mx,'string',cp(1,1));
%         set(handles.my,'string',cp(1,2));
%         set(handles.mvalue,'string',handles.cur_data(round(cp(1,1)),round(cp(1,2))));
%         %         handles.line_index=handles.line_index+1;
%         %         handles.line(handles.line_index,1)=cp(1,1);
%         %         handles.line(handles.line_index,2)=cp(1,2);
%     else
%         set(handles.figure1,'pointer','arrow');
%     end
% end

function NormalMotionFcn(hObject, eventdata, handles)
% if ~nargin
handles=guidata(gcbo);
pos=get(handles.OrgImg,'Position');
pt=get(handles.figure1,'CurrentPoint');
% pt=get(gcf,'CurrentPoint');
x=pt(1,1);
y=pt(1,2);

if (x>pos(1) & x<pos(1)+pos(3) & y>pos(2) & y<pos(2)+pos(4) )
    cp=get(handles.OrgImg,'CurrentPoint');

    if round(cp(1,2))> 0 & round(cp(1,2))<handles.frame_height & round(cp(1,1))>0 & round(cp(1,1))<handles.frame_width
        set(handles.figure1,'pointer','cross');
        set(handles.mx,'string',cp(1,1));
        set(handles.my,'string',cp(1,2));
        set(handles.mvalue,'string',handles.cur_data(round(cp(1,2)),round(cp(1,1))));
    else
        set(handles.figure1,'pointer','arrow');
    end

else
    set(handles.figure1,'pointer','arrow');
end



% --- Executes on mouse press over axes background.
function OrgImg_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to OrgImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function SlideWindow_Callback(hObject, eventdata, handles)
% hObject    handle to SlideWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SlideWindow as text
%        str2double(get(hObject,'String')) returns contents of SlideWindow as a double


% --- Executes during object creation, after setting all properties.
function SlideWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlideWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to SlideWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SlideWindow as text
%        str2double(get(hObject,'String')) returns contents of SlideWindow as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlideWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkstLASCA.
function checkstLASCA_Callback(hObject, eventdata, handles)
handles.enablestLASCA=~handles.enablestLASCA;
if handles.enablestLASCA==1
    set(handles.stLASCA,'enable','On');
else
    set(handles.stLASCA,'enable','off');
end
guidata(gcbo,handles);
% hObject    handle to checkstLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkstLASCA


% --- Executes on button press in stLASCA.
function stLASCA_Callback(hObject, eventdata, handles)
% hObject    handle to stLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data = handles.GroupstVeloImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['spatial & remporal speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;


% --- Executes on button press in checksLASCA.
function checksLASCA_Callback(hObject, eventdata, handles)
handles.enablesLASCA=~handles.enablesLASCA;
if handles.enablesLASCA==1
    set(handles.sLASCA,'enable','On');
else
    set(handles.sLASCA,'enable','off');
end
guidata(gcbo,handles);
% hObject    handle to checksLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checksLASCA


% --- Executes on button press in sLASCA.
function sLASCA_Callback(hObject, eventdata, handles)
% hObject    handle to sLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data = handles.GroupsImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['time-averaging spatial speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;


% --- Executes on button press in checktLASCA.
function checktLASCA_Callback(hObject, eventdata, handles)
handles.enabletLASCA=~handles.enabletLASCA;
if handles.enabletLASCA==1
    set(handles.tLASCA,'enable','On');
else
    set(handles.tLASCA,'enable','off');
end
guidata(gcbo,handles);
% hObject    handle to checktLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checktLASCA


% --- Executes on button press in tLASCA.
function tLASCA_Callback(hObject, eventdata, handles)
% hObject    handle to tLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data = handles.GrouptVeloImg(:,:,handles.currentgroup);
c=str2num(get(handles.SlideWindow,'string'));
[m,n]=size(handles.cur_data);
a=0;
for i=((c+1)/2):(m-(c-1)/2)
    for j=((c+1)/2):(n-(c-1)/2)
        for k=1:2:(2*c-1)
            for p=1:2:(2*c-1)
                a=a+handles.cur_data(i-(c-k)/2,j-(c-p)/2);
            end
        end
        handles.cur_data(i,j)=a/c^2;
        a=0;
    end
end
imagesc(handles.cur_data),title(['spatial-averaging temporal speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal


% --- Executes on button press in checknLASCA.
function checknLASCA_Callback(hObject, eventdata, handles)
% hObject    handle to checknLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.enablenLASCA=~handles.enablenLASCA;
if handles.enablenLASCA==1
    set(handles.nLASCA,'enable','On');
else
    set(handles.nLASCA,'enable','off');
end
guidata(gcbo,handles);
% Hint: get(hObject,'Value') returns toggle state of checknLASCA


% --- Executes on button press in nLASCA.
function nLASCA_Callback(hObject, eventdata, handles)
% hObject    handle to nLASCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data = handles.GroupnImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['normalize spatial speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;
% keyboard


% --- Executes on button press in checkaddpt.
function checkaddpt_Callback(hObject, eventdata, handles)
% hObject    handle to checkaddpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.enableaddpt=~handles.enableaddpt;
if handles.enableaddpt==1
    set(handles.addpt,'enable','On');
else
    set(handles.addpt,'enable','off');
end
guidata(gcbo,handles);
% Hint: get(hObject,'Value') returns toggle state of checkaddpt


% --- Executes on button press in addpt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 考虑到静态散射光
function addpt_Callback(hObject, eventdata, handles)
% hObject    handle to addpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data=handles.GroupaddVeloImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['temporal speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]);colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;

% --- Executes on button press in checkaddpk.
function checkaddpk_Callback(hObject, eventdata, handles)
% hObject    handle to checkaddpk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.enableaddpk=~handles.enableaddpk;
if handles.enableaddpk==1
    set(handles.addpk,'enable','On');
else
    set(handles.addpk,'enable','off');
end
guidata(gcbo,handles);
% Hint: get(hObject,'Value') returns toggle state of checkaddpk


% --- Executes on button press in addpk.
function addpk_Callback(hObject, eventdata, handles)
% hObject    handle to addpk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data = handles.GroupaddpkImg(:,:,handles.currentgroup);
imagesc(handles.cur_data),title(['spatial speckle contrast image: group ',...
    sprintf('%d',handles.currentgroup)]),colormap(gray);
handles.savefig=gcf;
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;


% --- Executes on button press in checkAVS.
function checkAVS_Callback(hObject, eventdata, handles)
% hObject    handle to checkAVS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.enableAVS=~handles.enableAVS;
if handles.enableAVS==1
    set(handles.AVS,'enable','On');
else
    set(handles.AVS,'enable','off');
end
guidata(gcbo,handles);
% Hint: get(hObject,'Value') returns toggle state of checkAVS


% --- Executes on button press in AVS.
function AVS_Callback(hObject, eventdata, handles)
% hObject    handle to AVS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cur_data = handles.GroupAVSVeloImg(:,:,handles.currentgroup);
imagesc(handles.cur_data);colormap(gray),title(['artery and vein separation image by contrast: group ',...
    sprintf('%d',handles.currentgroup)]);
handles.Cimage=getimage;

[m,n]=size(handles.cur_data);
frame_width=n;
frame_height=m;
handles.frame_width=frame_width;
handles.frame_height=frame_height;

guidata(gcbo,handles);
axis off;
axis equal;
