clear all
close all
clc
s=tf('s')
 G(1)=1/s/(s+1)^2;
 G(2)=1/s/(s^2-1);
 G(3)=s/(s^2+1);
 G(4)=1/s^2/(s-1);
 G(5)=s/(s^2-1);
 G(6)=(s-1)/(s^2*(s+10));
 for i=1:6
    s=tf('s')
    G(i)
    figure
    nyqlog(G(i))
end

