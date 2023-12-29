clear all, close all, clc
load 'C:\MATLABdata\allFaces.mat' %loading file
allPersons = zeros(n*6,m*6); %set variable allPersons to zeros matrix with #row*6 & #col*6, row&col in allFace.mat
count = 1; %set count = 1
for i=1:6 %Nested for loops
    for j=1:6 %Iterate row&col in allPersons matrix
        allPersons(1+(i-1)*n:i*n,1+(j-1)*m:j*m) ...
            = reshape(faces(:,1+sum(nfaces(1:count-1))),n,m); %update
        count = count + 1; %update count
    end
end %exit loops
figure(1), axes('position',[0  0  1  1]), axis off %print figure
imagesc(allPersons), colormap gray %set the color to be gray
%%
for person = 1:length(nfaces) %set person = 1, then for person to length of row of faces (dataset that built in matlab)
    subset = faces(:,1+sum(nfaces(1:person-1)):sum(nfaces(1:person)));  %set subset of faces
    allFaces = zeros(n*8,m*8); %update allFaces
    
    count = 1; %set the count = 1
    for i=1:8 %Nested for loops
        for j=1:8 %Iterate row&col in faces matrix
            if(count<=nfaces(person)) %if count <= row of faces of person
                allFaces(1+(i-1)*n:i*n,1+(j-1)*m:j*m) ... %then update allFaces
                    = reshape(subset(:,count),n,m);
                count = count + 1; %update count
            end
        end
    end 
    
    imagesc(allFaces), colormap gray    %print allfaces
end %exit loops