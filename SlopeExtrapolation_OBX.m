%% programatically extrapolate points to below a user specified elevation threshold. 
% 1/25/2019 by k.ames

% 3/31/2020 by J. Over
% Programatically create profiles from raw transect data and find transects
% that do not go below specified elevation threshold. Extrapolates points
% based on a slope threshold.

%----------------Inputs needed in Matlab_variables.xlsx file------------------------------------------------------------------------------------------------------
% A) SiteAbbv - Two letter site abbreviation, must match filename
% B) Path - path to where data files are located
% C) Export Path - where .shp file and .fig files will be output
% D) Start - (Landward or Seaward)- where transects start
% E) Sort Direction (Easting or Northing) site dependent, direction of transect sorting 
% G) ExportSHP - (Yes or No) - export a shape file of dune crest picks.
%                              WARNING this will overwrite any existing file. 
% I) Elevation Threshold - (number) - user determined elevation threshold.
%                           If no elevation on the transect crosses this threchold a point will be
%                           extrapolated and added into shapefile
%----------------Data Input---------------------------------------------------------------------------------------------------------------------------------------
% 1) XX_3D_transectdata.csv:(XX is site abbreviation) a csv with the 
%                           following columns: FID, Point, POINT_X, POINT_Y, elevation, Transect, INFO, Name
%                           note: INFO is for infilled points, Name is the survey file name ie 20141111
% 2) azimuth.csv: a csv with the X,Y, and azimuth of survey transects
%                 with the following columns: Shape_Length, OBJECTID, Transect, Start_X, Start_Y, End_X, End_Y, Azimuth
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
%%
close all
clear all
%-------------------USER EDITS BEGIN------------------------------------------------------------------------------------------------------------------------------
Site='CB'; %2 letter site abbreviation
slope_thresh=-.01; %change to control what is the minimum slope used to extrapolate
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
% load in Matlab variable sheet and find appropriate row based on site
% abbreviation
[~,UsrVar,~]=xlsread('Matlab_variables.xlsx');%,'A2:G50');
[~,~,Usr_dbthresh]=xlsread('Matlab_variables.xlsx','H2:I50');

[~, Usr_row]=max(strcmpi(Site,UsrVar(:,1)));
Usr_SiteVar=UsrVar(Usr_row,:);
Usr_elev_thresh= 0.4 %cell2mat(Usr_dbthresh(Usr_row,2));          %Fix
Usr_dbthresh= 2.5; %cell2mat(Usr_dbthresh(Usr_row,1));

%clear(UsrVar)
path=cell2mat(Usr_SiteVar(2));
path_exp=cell2mat(Usr_SiteVar(3));
Usr_start=cell2mat(Usr_SiteVar(4));
Usr_sort_dir=cell2mat(Usr_SiteVar(5));
Usr_shp=cell2mat(Usr_SiteVar(6)); %6 or 7?

path_filename=[path Site '_3D_transectdata.csv']; %name of excel file with survey data in it
Usr_azimuthfile=[path 'azimuth.csv']; %name of file with transect points in it start x, end x, start y, end y of transects

disp('Loading raw data from excel file...');

%adding first point on transect to transect point data, needed for consistant distance plotting
azimuthfile=csvread(Usr_azimuthfile,1,1);
azimuthfile(:,8)=repmat(.01,1,length(azimuthfile));
azimuth=azimuthfile(:,7); %There is a blank column in 7?

%reading survey data
raw=csvread(path_filename,1,1);
raw=sortrows(raw,[7 5 2]); 
Transect=raw(:,5); Transect2=unique(Transect).';
Sdate=raw(:,7);

dates=unique(Sdate).' %get unique survey dates from "Name" colume
timevar=num2str(dates(1,:));
seasons = numel(dates); %Counts the number of seasons/surveys
T_length=length(Transect);
num_transects=max(Transect);

dbthresh_z=[Usr_dbthresh,Usr_dbthresh];

TF_shp=strcmpi(Usr_shp,'yes');
TF_start=strcmpi(Usr_start,'seaward');
TF_sort = strcmpi(Usr_sort_dir,'Easting');
if TF_start==1 %depending on inland or seaward transect start get azimuth and start of transect
    azimuth=deg2rad(azimuth-180);
    tran_Xstart=azimuthfile(:,5);
    tran_Ystart=azimuthfile(:,6);
    tran_Xend=azimuthfile(:,3);
    tran_Yend=azimuthfile(:,4);
    if TF_sort==1
        sort_dir=3
    else 
        sort_dir=4
    end
else
    azimuth=deg2rad(azimuth-180);
    tran_Xstart=azimuthfile(:,3);
    tran_Ystart=azimuthfile(:,4);
    tran_Xend=azimuthfile(:,5);
    tran_Yend=azimuthfile(:,6);
    disp('Flipping azimuth')
    if TF_sort==1
        sort_dir=-3
    else 
        sort_dir=-4
    end
end
%% -------------------------------------------------split by season and transect and pick dune crest--------------------------------------------------------------
disp('Splitting data...'); %split by date
deletept = [];
for loop_i=dates %organize by dates
    index = find(Sdate==loop_i); %finds which index values are associated with the date
    genvarname('S_',  num2str(loop_i));  %create variable names with dates
    eval(['S_' num2str(loop_i) '=horzcat(raw(index,1:7),zeros(length(index),1));']) %splits raw data associated with each season into its own variable
    all_data(1,1:8)=zeros; %reset the all_data matrix
    all_data_orig(1,1:8)=zeros; %reset the all_data matrix
    all_ind=zeros(length(raw(:,2)),10);
    Transect_profile = [0];
    ext_point=zeros(1,8);
    for loop_jt=Transect2 %loop for each transect
        q= ['Transect = ', num2str(loop_jt) ,' season = ', num2str(loop_i) ,'']; %which transect and season is looping
        eval(['m1=S_', num2str(loop_i) ,'(:,5);'])
        index_trans=find(m1==loop_jt); %finds which index values within the season are associated with individual transects    

        if sum(index_trans==0)||(isempty(index_trans)==1)
            disp(['no data for ' q]) %display if missing data or transect
        else
            eval(['Transect_profile=S_', num2str(loop_i) ,'(index_trans,1:7);']);

            %set up and add fake elevation for beginning transect point
            fake_elev=repmat(1,1,1);
            azi_ind=find(azimuthfile(:,2)==loop_jt); 
            addon=[azimuthfile(azi_ind,1),tran_Xstart(azi_ind),tran_Ystart(azi_ind),fake_elev,azimuthfile(azi_ind,2),loop_i,loop_i];%start is landward 
            Transect_profile(end,:)=addon;
%----------------Projecting all points onto transect line--------------------------------------------------------------------------------------------------------------
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

%% -----------------sort and calculate distance between each point--------------------------------------------------------------------------------------------------------------
            Transect_profile=sortrows(Transect_profile,sort_dir); %sorts the data for a single transect by x or y
            for loop_ii=2:1:length(nonzeros(Transect_profile(:,2))) %calculate distance between points for all transects
                loop_tt=loop_ii-1;
                d(loop_ii,loop_jt)=sqrt((Transect_profile(loop_ii,2)-Transect_profile(loop_tt,2))^2 + (Transect_profile(loop_ii,3)-Transect_profile(loop_tt,3))^2);
                d(loop_ii,loop_jt)=d(loop_ii,loop_jt)+d(loop_tt,loop_jt);
            end   
            Transect_profile = horzcat(Transect_profile, d(1:loop_ii,loop_jt)); %concatenates the sorted data with its distance values 
            
            if sum(all_data_orig)==0
                all_data_orig = Transect_profile;
            else
                if all_data_orig(2,7)==loop_i
                    all_data_orig = vertcat(all_data_orig,Transect_profile);
                else
                    all_data_orig(1,1:8)=zeros;
                    all_data_orig = Transect_profile;
                end
            end   
%% -----------------Extrapolate point if no point on the profile crosses zero.---------------------------------------------------------------------------------------------------    
            elev_test=find(Transect_profile(:,4)<=Usr_elev_thresh);             %Leave for now and reintroduce threshold as ztarget below         
            if isempty(elev_test)==1
                %clear d1 d2 z1 z2 m m_neg m_neg2 
                
                ztest = Transect_profile(:,4);
                xtest = Transect_profile(:,8);
                ztarget = 0.4;                                                 %New user threshold
                %disp(['Threshold is ' num2str(ztarget)]);
                
                p = polyfit(xtest(end-4:end),ztest(end-4:end),1);              %fit straight line to last five points: output is slope and x-intercept
                xt = (ztarget-p(2))/p(1);                                      % rearrange line formula z = p(2)+p(1)*(x) to find x location of ztarget
                
                %Conditions for not extrapolating
                allBelow = find(ztest >= ztarget);                             %Checks to see if any points are above threshold, if not, don't put in point
                badExtra = find(xt < Transect_profile(end,8));                 %Checks to see if plotted point is behind last point, which would indicate profile already below threshold
                
                if p(1) >=0 %|| isempty(allBelow)==1                            %Or statement, p(1)>0 checks for no positive slopes
                    temp_ext_point=[1,NaN,NaN,ztarget,loop_jt,0,loop_i,NaN];   %Put in NaNs if no good data to extrapolate with
                    Transect_profile(end+1,:)=temp_ext_point;
                    Transect_profile(end+2:end,:)=[];
                    fprintf(' %.3f for profile %d during %d not plotted.\n', ztarget, loop_jt, loop_i);
                %elseif isempty(badExtra)==0
                    
                    
                    %temp_ext_point=[1,NaN,NaN,ztarget,loop_jt,0,loop_i,NaN];   %Put in NaNs if no good data to extrapolate with
                    %Transect_profile(end+1,:)=temp_ext_point;
                    %Transect_profile(end+2:end,:)=[];
                    %fprintf('profile %d during %d all below target.\n', loop_jt, loop_i);
                else
                    %make arrays to plot
                    xhat = [xtest(end-4):xtest(end)+20];
                    zhat = polyval(p,xhat);
 
                    %%find the slope in between points on the profile
                    z1=Transect_profile(1:end-1,4);
                    z2=Transect_profile(2:end,4);

                    d1=d(1:length(z1),loop_jt);
                    d2=d(2:length(z1)+1,loop_jt);
                    
                    %%This is old slope calculations for just 2 points
                    %m=(z2-z1)./(d2-d1); %calculate slope
                    %m(isnan(m))=0; %get rid of NaNs in calculation, replace with zeros
                    %ind=m<0; %find where slope is negative 
                    %m_neg=ind.*m; % output only negative slopes into a new matrix

                    %%find slopes below slope threshold (avoid using very minimum slopes) 
                    %i_thresh=m_neg<slope_thresh;
                    %m_neg2=m_neg.*i_thresh;
                    %[rowIdx,colIdx] = find(m_neg2); %get initial index of values
                    %v = accumarray(colIdx,rowIdx,[],@max).'; %get the max index (the last negative slope)
                
                    fprintf('Extrapolating slope below %.3f for profile %d during %d.\n', ztarget, loop_jt, loop_i);

                    %b=z1(v)-(m_neg2(v).*d1(v));
                    %%find x for given elevation threshold
                    %d_ext=(Usr_elev_thresh-b)./m_neg2(v);
                    %find last point before extrapolated point
                    d_diff=xt-(Transect_profile(end,8));  %%REPLACED d_ext with xt and d2(v) with distance end point
                    delta_X = d_diff.*sin(azimuth(loop_jt));
                    delta_Y = d_diff.*cos(azimuth(loop_jt));

                    New_X = Transect_profile(end,2) + delta_X;
                    New_Y = Transect_profile(end,3) + delta_Y;

                    temp_ext_point=[1,New_X,New_Y,ztarget,loop_jt,0,loop_i,xt];     %Puts the coordinates of the extrapolated point into a new matrix
                    Transect_profile(end+1,:)=temp_ext_point;
                    Transect_profile(end+2:end,:)=[];
                    
                    if sum(ext_point)==0
                        ext_point=temp_ext_point;
                    else
                        ext_point = vertcat(ext_point,temp_ext_point);
                    end  %end Dune
                end %end positive/negative if statment for good slope
%             else
%                 %attempt to find where crosses 0.4
%                 firstIndex = find(ztest <= 0.4,1, 'first');
%                 elevIndex = ztarget(firstIndex);
%                 distCross = xtest(firstIndex);
%                 if distCross < 50 | isempty(firstIndex)==1
%                     temp_ext_point=[1,NaN,NaN,ztarget,loop_jt,0,loop_i,NaN];   %Put in NaNs if no good data to extrapolate with
%                     Transect_profile(end+1,:)=temp_ext_point;
%                     Transect_profile(end+2:end,:)=[];
%                     fprintf('You will need to check profile %d during %d.\n',loop_jt, loop_i);
%                 else
%                     d_diff=distCross-(Transect_profile(end,8));  %
%                     delta_X = d_diff.*sin(azimuth(loop_jt));
%                     delta_Y = d_diff.*cos(azimuth(loop_jt));
% 
%                     New_X = Transect_profile(end,2) + delta_X;
%                     New_Y = Transect_profile(end,3) + delta_Y;
% 
%                     temp_ext_point=[1,New_X,New_Y,ztarget,loop_jt,0,loop_i,distCross];     %Puts the coordinates of the extrapolated point into a new matrix
%                     Transect_profile(end+1,:)=temp_ext_point;
%                     Transect_profile(end+2:end,:)=[];
%                     fprintf('Profile %d during %d crossed threshold at %d.\n',loop_jt, loop_i, distCross);
%                 end    
            end %end extrapolation if statement

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
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
        end %end check if transect is empty
    end %end JT for-loop  
    eval(['Ext_point_' num2str(loop_i) '=ext_point;']);
    TF_empty=sum(eval(['sum(Ext_point_' num2str(loop_i) ')']))==0;
    eval(['S_ext' num2str(loop_i) '=all_data;']); 
    eval(['S_' num2str(loop_i) '=all_data_orig;']);
%% ---export to feature class if set to yes-------------------------------------------------------------------------------------------------------------------------- 
%     if TF_shp==1 && TF_empty == 0
%         Geometry=repmat('Point',size(ext_point,1),1);
%         Ext_shp = struct('Geometry',cell(size(ext_point,1),1),'ID', cell(size(ext_point,1),1),'X', cell(size(ext_point,1),1),...
%         'Y', cell(size(ext_point,1),1), 'Z', cell(size(ext_point,1),1),'Transect',cell(size(ext_point,1),1),'INFO',cell(size(ext_point,1),1));       
%         for i1=1:size(ext_point,1) %set up structure
%             Ext_shp(i1).Geometry=Geometry(i1,:);
%             Ext_shp(i1).ID=ext_point(i1,1);
%             Ext_shp(i1).X=ext_point(i1,2);
%             Ext_shp(i1).Y=ext_point(i1,3);
%             Ext_shp(i1).Z=ext_point(i1,4);
%             Ext_shp(i1).Transect=ext_point(i1,5).';
%             Ext_shp(i1).INFO=ext_point(i1,6);
%         end
%         timevar=num2str(loop_i);
%         %dbfspec=makedbfspec(Ext_shp);
%         %filename3=[path Site '_extpoint_' timevar];   
%         %shapewrite(Ext_shp, filename3, 'DbfSpec', dbfspec)
%         
%         %disp(['Exporting extrapolated points for ' timevar]);
%     end
end %end date for loop

%% -------plotting data-------------------------------------------------------------------------------------------------------------------------- 
%                 plot(xtest,ztest);
%                 hold on
%                 plot(xhat,zhat,'--');
%                 plot(xt,.4,'ok')
% for loop_jt2=Transect2
%     f=figure(loop_jt2); %set up figure to call by transect and plot background stuff. 
%     x_dist=round(max(max(d1+50)),2,'significant'); 
%     y_dist=round(max(max(z1+1)),1,'significant');
%     %plot([0,x_dist],dbthresh_z,'g--')
%     hold on
%     plot([0,x_dist],[ztarget,ztarget],'b--')
%     %xlim([0,x_dist])
%     ylim([-1,y_dist])
%     title(strcat('Profile-',num2str(loop_jt2)))
%     xlabel('Distance(m)')
%     ylabel('Elevation(m)') 
% 
%     for loop_i2=dates
%         date1=num2str(loop_i2);
%         [d_s,loop_i]=sort(dates(1),8);
%         %eval(['t=Transect_',num2str(dates(1)),'(loop_i);']);
%         %figure(jt2); %recall transect figure
%         eval(['t_ind=find(S_',date1 '(:,5)==loop_jt2);']) %find correct transect to plot ## fixed so S_ doesn't have to be same size or in transect order. 
%         if isempty(t_ind)==0
%             a=eval(['plot(S_', date1 ,'(t_ind,8),S_', date1,'(t_ind,4));']);             
%             a.DisplayName=date1;
%             a.LineWidth=3;
%         end        
%         eval(['e_ind=find(S_ext',date1 '(:,5)==loop_jt2);'])   
%         eval(['d_ind=find(Ext_point_',date1 '(:,5)==loop_jt2);']) %find correct extrapolated point to plot
%         if isempty(d_ind)==0
%             h=eval(['plot(Ext_point_', date1 ,'(d_ind,8),Ext_point_', date1,'(d_ind,4));']); 
%             h.Color='black'; h.MarkerSize=15; h.Marker='.';h.DisplayName=strcat(date1,' Extrapolated'); h.LineStyle='none';   
%             e=eval(['plot(S_ext', date1 ,'(e_ind,8),S_ext', date1,'(e_ind,4));']);             
%             e.DisplayName=strcat(date1,' Extrapolated');
%             e.LineWidth=1;
%             e.LineStyle='--';
%             e.Color='black';
%         end
%     end
% end 

clear ans azi_ind b colIdx deletept Sdate Transect TF_shp TF_start timevar Usr_shp Usr_sort_dir Usr_start UsrVAR Usr_row
disp('End.');  
%------------------------------------------------------------------------------------------------------------------------------------------------------