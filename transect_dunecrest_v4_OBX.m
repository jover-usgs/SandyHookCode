%% programatically create profiles
% 7/14/2017 k.ames and k.butler
% edited 1/25/2019 by k.ames
% Programatically create profiles from raw transect data and choose maximum
% elevation on each profile. Plot all seasons on one figure for consistant
% dune crest id, custom data cursor (labeldtips4.m) for easy editing and
% picking. labeldtips4.m must be in same folder and in export matlab figure
% folder. 
% 1/18/2019 added projection to line for consistant plotting. 

%----------------Inputs needed in Matlab_variables.xlsx file------------------------------------------------------------------------------------------------------
% A) SiteAbbv - Two letter site abbreviation, must match filename
% B) Path - path to where data files are located
% C) Export Path - where .shp file and .fig files will be output
% D) Start - (Landward or Seaward)- where transects start
% E) UseComp? - (Yes or No) - use a .shp file to constrain where dune
%                             crests picks are made
% F) ExportSHP - (Yes or No) - export a shape file of dune crest picks.
%                              WARNING this will overwrite any existing file. 
% G) DuneThreshold - (number) - user determined elevation 
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------

%----------------Data Input---------------------------------------------------------------------------------------------------------------------------------------
% 1) XX_3D_transectdata.csv:(XX is site abbreviation) a csv with the 
%                           following columns: FID, Point, POINT_X, POINT_Y, elevation, Transect, INFO, Name
%                           note: INFO is for infilled points, Name is the survey file name ie 20141111
% 2) azimuth.csv: a csv with the X,Y, and azimuth of survey transects
%                 with the following columns: Shape_Length, OBJECTID, Transect, End_X, End_Y, Start_X, Start_Y, Azimuth
% 3) compartments.shp: a shapefile containing the verticies of the compartments used with the following attributes: 
%                      Point_X, Point_Y, Order
%                      note: the order field is needed for the compartment polygon to draw
%                      correctly. The seaward line of points order =1, landward side
%                      order=2. Then the verticies are sorted by these field then by
%                      northing and easting. To change sort order go to line 87 & 90. the
%                      command is "sortrows(shp_xy2,1);"
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
% IMPORTANT! currently written for 10 season of data, if greater than 10 seasons must manually change  labeldtips4.m to read correct number of seasons as input.
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%%
close all
clear all
%-------------------USER EDITS BEGIN------------------------------------------------------------------------------------------------------------------------------
Site='CB'; %2 letter site abbreviation
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
% load in Matlab variable sheet and find appropriate row based on site
% abbreviation
[~,UsrVar,~]=xlsread('Matlab_variables.xlsx'); %,'A2:G50');
[~,~,Usr_dbthresh]=xlsread('Matlab_variables.xlsx','H2:I50');

[~, Usr_row]=max(strcmpi(Site,UsrVar(:,1)));
Usr_SiteVar=UsrVar(Usr_row,:);
Usr_elev_thresh= 0 %cell2mat(Usr_dbthresh(Usr_row,2));
Usr_dbthresh=2.5 %cell2mat(Usr_dbthresh(Usr_row,1));

%clear(UsrVar) - the following can be changed
path=cell2mat(Usr_SiteVar(2));
path_exp=cell2mat(Usr_SiteVar(3));
Usr_start=cell2mat(Usr_SiteVar(4));
Usr_sort_dir=cell2mat(Usr_SiteVar(5));
comp_use=cell2mat(Usr_SiteVar(8));
Usr_shp=cell2mat(Usr_SiteVar(6));

elev_thresh=-0.335; %elevation of base raster in /Infill/ folder
NAVD=0; %NAVD88 threshold
filename=[path Site '_3D_transectdata.csv']; %name of excel file with survey data in it
Usr_azimuthfile=[path 'azimuth.csv']; %name of file with transect points in it start x, end x, start y, end y of transects
compfilename=[path 'CompartmentPoints.shp']; %name of shape file with compartment bounding points in it (feature verticies to points in Arcmap) 
 
compTF=strcmp(comp_use,'yes');
if compTF==1
    compfilename=[path 'CompartmentPoints.shp']; 
end
disp('Loading raw data from excel file...');

%adding first point on transect to transect point data, needed for consistant distance plotting
azimuthfile=csvread(Usr_azimuthfile,1,1);
azimuthfile(:,8)=repmat(.01,1,length(azimuthfile));
azimuth=azimuthfile(:,7);

%reading survey data
raw=csvread(filename,1,1);
raw=sortrows(raw,[7 5 2]); Point=raw(:,1);
x=raw(:,2); y=raw(:,3); z=raw(:,4);
Transect=raw(:,5); Transect2=unique(Transect).';
Rdate=raw(:,6); Sdate=raw(:,7);

dates=unique(Sdate).' %get unique survey dates from "Name" colume
timevar=num2str(dates(1,:));
seasons = numel(dates); %Counts the number of seasons/surveys
T_length=length(Transect);
num_transects=max(Transect);

dbthresh_z=[Usr_dbthresh,Usr_dbthresh];
NAVD_z=[NAVD,NAVD];
if compTF==1
    %loading in compartment verticies
    disp('Loading data from shapefile...');
    dune_limit = shaperead(compfilename); 
    %extract x and y values from the shapefile
    shp_xy(:,1)=extractfield(dune_limit,'X').'; 
    shp_xy(:,2)=extractfield(dune_limit,'Y').'; 
    %this is a field in the shapefile seasward is 1 landward is 2. needed to
    %sort the polygon, can change the order of the sort in lines 86 and 88
    shp_xy(:,3)=extractfield(dune_limit,'Order_').';

    shp_xy=sortrows(shp_xy,3); %sort by order field
    in_order=find(shp_xy(:,3)==1); %find the seaward points

    shp_xy1=shp_xy(1:max(in_order),:); %split into seaward points
    shp_xy1=sortrows(shp_xy1,-1); %and sort by decending northing (change to 2 if sorting by easting)

    shp_xy2=shp_xy(max(in_order)+1:end,:); %split into landward points
    shp_xy2=sortrows(shp_xy2,1); %and sort by ascending northing to complete the shape
    shp_xy=[shp_xy1;shp_xy2]; %put data back together
    figure
    fill(shp_xy(:,1),shp_xy(:,2),'r');
end

TF_shp=strcmpi(Usr_shp,'yes');
TF_start=strcmpi(Usr_start,'landward');

if TF_start==1 %depending on inland or seaward transect start get azimuth and start of transect
    azimuth=deg2rad(azimuth);
    tran_Xstart=azimuthfile(:,5);
    tran_Ystart=azimuthfile(:,6);
    tran_Xend=azimuthfile(:,3);
    tran_Yend=azimuthfile(:,4);
else
    azimuth=deg2rad(azimuth-180);
    tran_Xstart=azimuthfile(:,4);
    tran_Ystart=azimuthfile(:,5);
    tran_Xend=azimuthfile(:,2);
    tran_Yend=azimuthfile(:,3);
    disp('Flipping azimuth')
end
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%if you need to check shape of compartment file uncomment two lines below
%figure
%fill(shp_xy(:,1),shp_xy(:,2),'r');
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------split by season and transect and pick dune crest--------------------------------------------------------------
disp('Splitting data...'); %split by date
deletept = [];
for loop_i=dates %organize by dates
    index = find(Sdate==loop_i); %finds which index values are associated with the date
    genvarname('S_',  num2str(loop_i));  %create variable names with dates
    eval(['S_' num2str(loop_i) '=raw(index,1:7);']);  %splits raw data associated with each season into its own variable
    all_data(1,1:8)=zeros; %reset the all_data matrix
    all_ind=zeros(length(x),10);
    Transect_profile = [0];
    DuneCrest=zeros(1,8);
        for loop_jt=Transect2 %loop for each transect
            q= ['Transect = ', num2str(loop_jt) ,' season = ', num2str(loop_i) ,'']; %which transect and season is looping
            eval(['m1=S_', num2str(loop_i) ,'(:,5);'])
            index_trans=find(m1==loop_jt); %finds which index values within the season are associated with individual transects    
            %temp_ind=index_trans;

            if sum(index_trans==0)||(isempty(index_trans)==1)
                disp(['no data for ' q]) %display if missing data or transect
            else
                if loop_jt==1
                   all_ind(1:length(index_trans),1)=index_trans;
                else
                    all_ind(1:length(index_trans),loop_jt) = index_trans;
                end 

                eval(['Transect_' num2str(loop_i) '=all_ind;']); %make an index variable for each season
                eval(['Transect_profile=S_', num2str(loop_i) ,'(index_trans,1:7);']);
                
                %set up and add fake elevation for beginning transect point
                %for consistant plotting. 
                fake_elev=repmat(.01,1,1);
                azi_ind=find(azimuthfile(:,2)==loop_jt); %Was previously 2, may change if figure out correct order
                if TF_start==1
                    addon=[azimuthfile(azi_ind,1),azimuthfile(azi_ind,3),azimuthfile(azi_ind,4),fake_elev,azimuthfile(azi_ind,2),loop_i,loop_i];%start is seaward
                else
                    addon=[azimuthfile(azi_ind,1),azimuthfile(azi_ind,5),azimuthfile(azi_ind,6),fake_elev,azimuthfile(azi_ind,2),loop_i,loop_i];%start is landward
                end
                Transect_profile(end,:)=addon;
                
%----------------Projecting all points onto transect line---------------------------------------
                delta_Xp = sin(azimuth(loop_jt));
                delta_Yp = cos(azimuth(loop_jt));

                %difference of point with transect origin
                dx = Transect_profile(:,2)-tran_Xstart(loop_jt);
                dy = Transect_profile(:,3)-tran_Ystart(loop_jt);

                %Position of projection on line, using dot product
                xs = (((dx.*delta_Xp)+(dy.*delta_Yp))/((delta_Xp^2)+(delta_Yp^2)))*delta_Xp;
                ys = (((dx.*delta_Xp)+(dy.*delta_Yp))/((delta_Xp^2)+(delta_Yp^2)))*delta_Yp;

                %convert point to position
                proj_x=tran_Xstart(loop_jt)+xs;
                proj_y=tran_Ystart(loop_jt)+ys;

                Transect_profile(:,2)=proj_x;
                Transect_profile(:,3)=proj_y;

%-----------------End Projection----------------------------------------------------------------
                %sort and calculate distance between each point
                Transect_profile=sortrows(Transect_profile,-3); %sorts the data for a single transect by x or y
                for loop_ii=2:1:length(nonzeros(Transect_profile(:,2))) %calculate distance between points for all transects
                                loop_tt=loop_ii-1;
                                d(loop_ii,loop_jt)=sqrt((Transect_profile(loop_ii,2)-Transect_profile(loop_tt,2))^2 + (Transect_profile(loop_ii,3)-Transect_profile(loop_tt,3))^2);
                                d(loop_ii,loop_jt)=d(loop_ii,loop_jt)+d(loop_tt,loop_jt);
                end   

                Transect_profile = horzcat(Transect_profile, d(1:loop_ii,loop_jt)); %concatenates the sorted data with its distance values                
                
                    if sum(all_data)==0
                         all_data = Transect_profile;
                    else
                        if all_data(2,7)==loop_i
                            all_data = vertcat(all_data,Transect_profile);
                        else
                            all_data(1,1:8)=zeros;
                            all_data = Transect_profile;
                        end
                    end   

                if compTF==1
                    in = inpolygon(Transect_profile(:,2),Transect_profile(:,3),shp_xy(:,1),shp_xy(:,2)); %find only data within compartments. 
                    z_lim=Transect_profile(:,4).*in;
                    if max(z_lim>0) %if there is a point inside the compartment polygon then pick dune crest
                        duneindex=find(Transect_profile(:,4)==max(z_lim));
                        duneindex=duneindex(1);
                        Dune_max=Transect_profile(duneindex,:);
                        if Dune_max(:,4)>=Usr_dbthresh
                            temp_DuneCrest=Dune_max;
                         else
                             temp_dunecrest=zeros(1,8);
                        end
                    end
                else
                    [~,duneindex]=max(Transect_profile(:,4));
                    duneindex=duneindex(1);
                    Dune_max=Transect_profile(duneindex,:);
                    if Dune_max(:,4)>=Usr_dbthresh
                        temp_DuneCrest=Dune_max;
                    else
                        temp_DuneCrest=Dune_max;
                        temp_DuneCrest(:,4)=NaN;
                    end
                end
                
                if sum(DuneCrest)==0
                    DuneCrest=temp_DuneCrest;
                else
                    DuneCrest = vertcat(DuneCrest,temp_DuneCrest);
                end  %end Dune
            end            
         end %end JT for-loop  
       
    eval(['S_' num2str(loop_i) '=all_data;']);
    eval(['Dune_' num2str(loop_i) '=DuneCrest;']);
      
    if TF_shp==1         %%export to feature class if set to yes
        Geometry=repmat('Point',length(DuneCrest),1);
        Dune_shp = struct('Geometry',cell(length(DuneCrest),1),'ID', cell(length(DuneCrest),1),'X', cell(length(DuneCrest),1),...
            'Y', cell(length(DuneCrest),1), 'Z', cell(length(DuneCrest),1),'Transect',cell(length(DuneCrest),1),'INFO',cell(length(DuneCrest),1));       
        for i1=1:length(DuneCrest) %set up structure
                Dune_shp(i1).Geometry=Geometry(i1,:);
                Dune_shp(i1).ID=DuneCrest(i1,1);
                Dune_shp(i1).X=DuneCrest(i1,2);
                Dune_shp(i1).Y=DuneCrest(i1,3);
                Dune_shp(i1).Z=DuneCrest(i1,4);
                Dune_shp(i1).Transect=DuneCrest(i1,5).';
                Dune_shp(i1).INFO=DuneCrest(i1,6);
        end
                timevar=num2str(loop_i);
                %dbfspec=makedbfspec(Dune_shp);
                %filename3=[path Site '_Dune_' timevar];   
                %shapewrite(Dune_shp, filename3, 'DbfSpec', dbfspec) 
                %disp(['Exporting Dune crest points for ' timevar]);
    end
end 

%append zeros in case last transects missing. For now this is only needed
%for broken data tips function. need a better solution. 
c=zeros(length(dates),1);
for loop_i3=dates
    %find widest matrix
    %[r,temp_c]= size(S_date);
    date2=num2str(loop_i3);
    eval(['[~,temp_c]= size(Transect_', date2 ');'])
    if sum(c) ==0
        c=temp_c;
    else
        c=vertcat(c,temp_c);
    end
    max_c=max(c);
    
    if temp_c<max(c)
        eval(['Transect_',date2 '(:,temp_c+1:max_c)=0;'])
    end
end
%% -----------------------------------------------plot all transects with all 10 seasons and dune crest pics -----------------------------------------------------
%plotting with custom data cursor
for loop_jt2=Transect2
    f=figure(loop_jt2); %set up figure to call by transect and plot background stuff. 
%     set(f,'Position', [0 0 1100 500])
%     set(f, 'PaperUnits', 'inches');
%     set(f, 'PaperSize', [9 5]);
%     set(f, 'PaperPositionMode', 'manual');
%     set(f, 'PaperPosition', [0 0 9 5]);
    x_dist=round(max(max(d+50)),2,'significant'); %the fake point has been added and it throwing off the plot
    y_dist=round(max(max(z+1)),1,'significant');
    plot([0,x_dist],dbthresh_z,'g--')
    hold on
    plot([0,x_dist],NAVD_z,'b--')
    xlim([0,x_dist])
    ylim([-1,y_dist])
    title(strcat('Profile-',num2str(loop_jt2)))
    xlabel('Distance(m)')
    ylabel('Elevation(m)') 
    

   for loop_i2=dates
       date1=num2str(loop_i2);
       %setting up colors for plots to loop through nicely
        set(groot,'DefaultAxesColorOrder',[0.843137255 0.109803922 0.152941176; 0 0 0;  0.956862745 0.42745098 0.262745098; 0 0 0;...
        0.992156863 0.682352941 0.380392157; 0 0 0; 0.996078431	0.878431373	0.564705882; 0 0 0; 0	0.77254902	1;...
        0 0 0; 0.670588235 0.850980392 0.91372549; 0 0 0; 0.454901961 0.678431373	0.819607843; 0 0 0; 0.270588235 0.458823529 0.705882353; 0 0 0;...
        0.5 0.5 0.5; 0 0 0; .25 .25 0]);
%         set(groot,'DefaultAxesColorOrder',[0.843137255 0.109803922 0.152941176;...
%          0.992156863 0.682352941 0.380392157; 0	0.77254902	1;...
%          0.454901961 0.678431373	0.819607843; ...
%          0.5 0.5 0.5;  .25 .25 0])
        
        [d_s,loop_i]=sort(dates(1),8);
        eval(['t=Transect_',num2str(dates(1)),'(loop_i);']);
        %figure(jt2); %recall transect figure
        eval(['t_ind=find(S_',date1 '(:,5)==loop_jt2);']) %find correct transect to plot ## fixed so S_ doesn't have to be same size or in transect order. 
        if isempty(t_ind)==0
            a=eval(['plot(S_', date1 ,'(t_ind,8),S_', date1,'(t_ind,4));']);             
            a.DisplayName=date1;
            a.LineWidth=1.5;
            %legend('','',dates(1),dates(2),dates(3))%,dates(4),dates(5), dates(6))
        end
        eval(['d_ind=find(Dune_',date1 '(:,5)==loop_jt2);']) %find correct dune crest to plot
        if isempty(d_ind)==0
            h=eval(['plot(Dune_', date1 ,'(d_ind,8),Dune_', date1,'(d_ind,4));']); 
            h.Color='black'; h.MarkerSize=15; h.Marker='.';h.DisplayName=date1; h.LineStyle='none';
        end
        
   end
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------   
%% set up custom datacursor for easy arcmap use
% set up custom datacursor for easy arcmap use
% IMPORTANT! currently written for  up to 10 season of data, if greater than 10 seasons must manualy change labeldtips4.m to
% read correct number of seasons as input.
hdt = datacursormode;
datacursormode off
set(hdt,'DisplayStyle','window');

%testing data cursor
if seasons>=1
    season1 = eval(['S_', num2str(dates(1)),'(nonzeros(Transect_', num2str(dates(1)),'(:,loop_jt2)),:);']);
    season2=0; season3=0; season4=0; season5=0; season6=0; season7=0; season8=0; season9=0; season10=0;  
end
if seasons>=2
    season2 = eval(['S_', num2str(dates(2)),'(nonzeros(Transect_', num2str(dates(2)),'(:,loop_jt2)),:);']);
    season3=0; season4=0; season5=0; season6=0; season7=0; season8=0; season9=0; season10=0;
end
if seasons>=3
    season3 = eval(['S_', num2str(dates(3)),'(nonzeros(Transect_', num2str(dates(3)),'(:,loop_jt2)),:);']);
    season4=0; season5=0; season6=0; season7=0; season8=0; season9=0; season10=0;
end
if seasons>=4
    season4 = eval(['S_', num2str(dates(4)),'(nonzeros(Transect_', num2str(dates(4)),'(:,loop_jt2)),:);']);
    season5=0; season6=0; season7=0; season8=0; season9=0; season10=0;
end
if seasons>=5
    season5 = eval(['S_', num2str(dates(5)),'(nonzeros(Transect_', num2str(dates(5)),'(:,loop_jt2)),:);']);
    season6=0; season7=0; season8=0; season9=0; season10=0;
end
if seasons>=6
    season6 = eval(['S_', num2str(dates(6)),'(nonzeros(Transect_', num2str(dates(6)),'(:,loop_jt2)),:);']);
    season7=0; season8=0; season9=0; season10=0;
end
if seasons>=7
    season7 = eval(['S_', num2str(dates(7)),'(nonzeros(Transect_', num2str(dates(7)),'(:,loop_jt2)),:);']);
    season8=0; season9=0; season10=0;
end
if seasons>=8
    season8 = eval(['S_', num2str(dates(8)),'(nonzeros(Transect_', num2str(dates(8)),'(:,loop_jt2)),:);']);
    season9=0; season10=0;
end
if seasons>=9
    season9 = eval(['S_', num2str(dates(9)),'(nonzeros(Transect_', num2str(dates(9)),'(:,loop_jt2)),:);']);
    season10=0;
end
if seasons>=10
    season10 = eval(['S_', num2str(dates(10)),'(nonzeros(Transect_', num2str(dates(10)),'(:,loop_jt2)),:);']);
end

set(hdt,'UpdateFcn',{@labeldtips4,season1,season2,season3,season4,season5,season6,season7,season8,season9,season10})
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%% export figure
fig_title=(['Profile-', num2str(loop_jt2), '']);
figure_name=[path fig_title];
%%orient(gcf,'landscape')
saveas(gcf,figure_name); %saving as matlab fig
%%saveas(gcf,figure_name,'pdf'); %saving as matlab fig
%%print(figure_name,'-dpdf','-fillpage')
end
%[pointslist,xselect,yselect] = selectdata
disp('End.');  
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------