function [pos,n]=group(x,crit,excl)
%function [pos,n]=group(x,crit,excl)
%determine the number of groups of an increasing vector with as criterium that it's group 
%all the elements which are distant of striclty less than 'crit' elements. It does not consider
%group that have less than 'excl' elements.
%pos gives the positions related to the vector x of the begins and ends of each group
%example group([1 2 4 6 10 11 15 16 18],2,2) gives n=2 and pos=[1 4 7 9];

[r l]=size(x);
if r~=1
  x=x';
end

temp=x(2:end)-x(1:end-1);
if min(temp)<1|min([r l])>1 %not increasing vector | a single element
  pos=[];
  n=0;
else
  s=find(temp>crit);
  pos(1)=1;
  t=2;
  for i=1:length(s)
     if (s(i)-pos(t-1))>=excl
        pos(t)=s(i);
        t=t+1;
        pos(t)=s(i)+1;
        t=t+1;
     else
        pos=setdiff(pos,pos(t-1));
        pos(t-1)=s(i)+1;
     end
  end
  if (length(x)-pos(t-1))>=excl
     pos(t)=length(x);
  else
     pos=setdiff(pos,pos(t-1));
     t=t-1;
  end   
  n=length(pos)/2;
end
