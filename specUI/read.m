%figure;
for i=1:25
    filename=['T',num2str(i)];
    load(filename);
    imagesc(A);
    axis off;
    %colormap(gray);
    %caxis([0 400]); %%gray
    caxis([0.02 0.30]); %%color
    %s=0;
    %for i=1:5
        %h=improfile;
        %[m,n]=size(h);
        %s=s+m;
    %end    
    re=getrect;
end