%plot vecteur 3d 

hold off;
colordef black;

solvent1=load('Wat_newprog.out');
solvent2=load('Met_newprog.out');

if 1 == 1
  
%solvent 1
dens1 = solvent1(:,2);  %extrait the probability presence at this place
vx1 = solvent1(:,3) ;   %extrait des vecteurs
vy1 = solvent1(:,4) ;
vz1 = solvent1(:,5) ;
gx1 = solvent1(:,6) ;   %extrait de leur origine
gy1 = solvent1(:,7) ;
gz1 = solvent1(:,8) ;

%solvent 2
dens2 = solvent2(:,2);  %extrait the probability presence at this place
vx2 = solvent2(:,3) ;   %extrait des vecteurs
vy2 = solvent2(:,4) ;
vz2 = solvent2(:,5) ;
gx2 = solvent2(:,6) ;   %extrait de leur origine
gy2 = solvent2(:,7) ;
gz2 = solvent2(:,8) ;

%solvent 1 
gxx1 =  [zeros(size(gx1))+gx1];
gyy1 =  [zeros(size(gy1))+gy1];
gzz1 =  [zeros(size(gz1))+gz1];
xx1 = [gxx1 vx1] ;
yy1 = [gyy1 vy1] ;
zz1 = [gzz1 vz1] ;

%solvent 2
gxx2 =  [zeros(size(gx2))+gx2];
gyy2 =  [zeros(size(gy2))+gy2];
gzz2 =  [zeros(size(gz2))+gz2];
xx2 = [gxx2 vx2] ;
yy2 = [gyy2 vy2] ;
zz2 = [gzz2 vz2] ;

%criteria selection
crit_pol = round(10*get(HF9,'Value'))/100;
crit_dens = round(1*get(HF10,'Value'))/10000;

norm_speed1 = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
low_speed1 = find(norm_speed1<crit_pol);
low_dens1 = find(dens1<crit_dens);
low_dens2 = find(dens2<crit_dens);

vx1(low_speed1) = NaN;
vy1(low_speed1) = NaN;
vz1(low_speed1) = NaN;

gxx1(low_speed1) = NaN;
gyy1(low_speed1) = NaN;
gzz1(low_speed1) = NaN;

vx1(low_dens1) = NaN;
vy1(low_dens1) = NaN;
vz1(low_dens1) = NaN;

gxx1(low_dens1) = NaN;
gyy1(low_dens1) = NaN;
gzz1(low_dens1) = NaN;

norm_speed2 = sqrt(vx2.^2 + vy2.^2 + vz2.^2);
low_speed2 = find(norm_speed2<crit_pol);

vx2(low_speed2) = NaN;
vy2(low_speed2) = NaN;
vz2(low_speed2) = NaN;

gxx2(low_speed2) = NaN;
gyy2(low_speed2) = NaN;
gzz2(low_speed2) = NaN;

vx2(low_dens2) = NaN;
vy2(low_dens2) = NaN;
vz2(low_dens2) = NaN;

gxx2(low_dens2) = NaN;
gyy2(low_dens2) = NaN;
gzz2(low_dens2) = NaN;


    cla
    [x,y,z] = sphere(20);
    fid = fopen('solute_pos.dat','r');
    line = fgetl(fid);
    while isstr(line)
            switch(line(8))
            case  'H', color = [1.0 1.0 1.0]; r = 0.8; 
            case  'C', color = [0.0 1.0 0.0]; r = 1.2; 
            case  'O', color = [1.0 0.0 0.0]; r = 1.2; 
            case  'N', color = [1.0 0.3 1.0]; r = 1.0; 
            otherwise, color = [0.6 0.6 0.6]; r = 1.2; 
            end
            c =  sscanf(line(1:32),'%s %i %i %i');
            surf(x.*r+c(2),y.*r+c(3),z.*r+c(4),'FaceColor',color,...
                'FaceLighting','gouraud','EdgeColor','none',...
                'AmbientStrength',0.65);
            light('Position',[-1 -1 -2]);
            hold on;
            line = fgetl(fid);
    end
    fclose(fid);

  
      %plot3(gxx1,gyy1,gzz1,'b.') ;
      %plot3(gxx2,gyy2,gzz2,'y.') ;
      
      quiver3_perso(gxx1,gyy1,gzz1,vx1,vy1,vz1);
      hold on;
      quiver3(gxx2,gyy2,gzz2,vx2,vy2,vz2);
      grid on;
      
     
xlabel('x direction'),ylabel('y direction'),zlabel('z direction');
title('density presence and orientation of solvent');
      hold off;

end