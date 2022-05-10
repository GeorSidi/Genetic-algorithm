function ask3
pkg load symbolic;

e=10^-6;
Dmin = -1024;         ##a
Dmax = 1023;          ##b
I = Dmax - Dmin;
k=1+I/e;
n=0;
while(1)
   if(k > 2^(n-1) && k <= 2^n)
    break;
   endif
   n++;
endwhile

m=n+n;                ##the number of bit we need cause we have 2 inputs with n bits each
p=500;                ##the number of population

############initialization######################################
x=randi([0 1], p ,m);
x(2,:);
xnew=zeros(p,m);      ##array for population after roulette
x_start=x;

x1_best=0;            ##x1,x2 for the min value of f
x2_best=0;            ##
f_min=1000;    ##a random start 

xb=zeros(2,p);        ##Arrays for x,z and finest f values
zb=zeros(2,p);        ##using for the encoding and decoding process
ff=zeros(1,p);        ##

############starting the loop######################################
for loop=1:100        ##we can set here the number of generations(100)        

############encoding process to take 0 1 bit arrays for x1,x2 #####
for i=1:p
  for j=0:n-1
    zb(1,i)=zb(1,i)+x(i,j+1)*(2^(j));
    zb(2,i)=zb(2,i)+x(i,j+1+n-1)*(2^(j));
  endfor
xb(1,i)=Dmin+zb(1,i)*(I/( (2^n) - 1 ));
xb(2,i)=Dmin+zb(2,i)*(I/( (2^n) - 1 ));
ff(i)=xb(1,i)^2 + xb(2,i);  ##we are minimized the fx= x1^2 + x2

if(ff(i)<f_min)       ##we are keeping the less value of f
  f_min=ff(i);        ##for every time at start of the loop
  x1_best=xb(1,i);    ##the same for x1 
  x2_best=xb(2,i);    ##the same for x2
endif
endfor
############end of initialations################################

############Roulette Selection##################################
fworst=max(ff);
b=find(fworst==ff);                   ##finding the worst value of f

allf=sum(ff);
totalfitness = (p*fworst) - allf;

fitness_values=zeros(1,p);
selection_possibility=zeros(1,p);

r=rand(1,1);
chosenone=1;

for i=1:p
  fitness_values(1,i)=fworst-ff(1,i);
  selection_possibility(1,i)= fitness_values(1,i)/totalfitness ;
  
  k=sum(selection_possibility(1:i)); ##k helps us to find the place of
  if(k<r)                            ##the chance,to pick a parent for 
    chosenone++;                     ##process
  endif  
xnew(i,:)=x(chosenone,:);            ##fixing the array
endfor

########################end roylette###########################

##################cross########################################
xcross=xnew;                         ##array for population after crossing
cp=4/10;
lastforcross=1;
forcrossed=[];
ncross=0;

for i=1:p
  rcross=rand(1,1);                  ##random value for crossing
  if(rcross>cp)   
    forcrossed=[forcrossed i];
    ncross++;
  else                               ##taking 1 more parent to cross 
    lastforcross++;                  ##if the value of selected parents
  end                                ##is odd
endfor
if(mod(ncross,2)>0)
  forcrossed=[forcrossed lastforcross];
endif

crossing=length(forcrossed);        ##we are using 1 point croossing in the 
wh=1;                               ##middle of each crossed parents
while(wh<crossing)
  xcross(forcrossed(wh),:)= [xnew(forcrossed(wh),1:n) , xnew(forcrossed(wh+1),n+1:end)];
  xcross(forcrossed(wh+1),:)= [xnew(forcrossed(wh+1),1:n) , xnew(forcrossed(wh),n+1:end)];
  wh+=2;
end
##########end croos###########################################

########mutations#############################################
xmut=xcross;                        ##array for population after mutations
mp=10/175;                          ##chance of mutation
for i=1:p
for j=1:m
  mutr=rand(1,1);                   ##if a random value is bigger than 
  if(mutr<mp)                       ##mp the bit will change
  if(xmut(i,j)==1)                  ##(for 1 to 0 or 0 to 1)
  xmut(i,j)=0;
  else
  xmut(i,j)=1;
  end
  endif
endfor
endfor

#####################End of Mutation#############################
                                   ##the end of every loop cycle to
x=xmut;                            ##set the new better population  
endfor                             ##for the next lopp cycle
##################End of the loop################################

##################Printing the Results###########################                           
x1_best
x2_best
f_min

endfunction
