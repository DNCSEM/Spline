function yi=splineyy(x,y,xi)
% This is cubic spline function, created by YangYang
%yi is the output vector, while x,y,xi is the input vector.
    if size(x,1)==1
        x=x';
    end
    if size(y,1)==1
        y=y';
    end
    if size(xi,1)==1
        xi=xi';
    end
    yi=zeros(length(xi),1);
    n=length(x);
    index=1:length(x);
    x1=[x,index'];
    xi1=[xi,zeros(length(xi),1)];
    xSort=sortrows([x1;xi1],1);
    xSort(end-2:end)=sortrows(xSort(end-2:end),[1 -2]);
    [r_sort,c_sort]=find(xSort(:,2)~=0);
     
    matrix=toeplitz([4,1,sparse(zeros(1,n-2))]);
    matrix([1,end])=2;
    
    new_y=[y(3:end);y(1:2)];
    diff_y=new_y-y;
    vector=[y(2)-y(1);diff_y(1:end-2);y(end)-y(end-1)];
    
    D=3*vector\matrix; 
    begin=0;
    for i=1:n-1
        num=r_sort(i+1)-r_sort(i);
        if num==1
            continue;
        end
        r=begin+1:begin+num-1;
        a=y(i);
        b=D(i);
        c=3*(y(i+1)-y(i))-2*D(i)-D(i+1);
        d=2*(y(i)-y(i+1))+D(i)+D(i+1);
        xi_r=(xi(r)-x(i))/(x(i+1)-x(i));
        yi(r)=a+b*xi_r+c*xi_r.^2+d*xi_r.^3;
        begin=begin+num-1;
    end
    plot(x,y,'ro',xi,yi);grid
end