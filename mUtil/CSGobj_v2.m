classdef CSGobj_v2
    properties
        PL; %Polygon set
        LS; %line Set
        PS; %point Set
        tol; %tolerance
        %GRD; %grid object for faster search
        dummy; %empty variable for transfering around values
    end
    
    methods

        %initialize object
        function CSG=CSGobj_v2(Dim,Npoly,Nline,Npoints,usertol)
            %Dim dimension of points
            CSG.PL=PolySet_v2(Npoly,usertol);
            CSG.LS=LineSet_v2(Nline,usertol);
            CSG.PS=PointSet_v2(Npoints,Dim,usertol);
            CSG.tol=usertol;
        end
        
        
        %input methods
        function CSG=readshapefile(CSG,SH,dis,offset)
            shapetype=SH(1,1).Geometry;
            switch shapetype
                case 'Polygon'
                    if CSG.PL.Npl>0
                        %a polygon has already be defined 
                        SH(1,1).Geometry='Line';
                        CSG=CSG.readshapefile(SH);
                        return;
                    end
                    CSG.dummy=1;
                    for ii=1:size(SH,1)
                        [Xs, Ys]=msim_polysplit(SH(ii,1).X,SH(ii,1).Y);
                        for jj=1:size(Xs,1)
                            if jj>size(CSG.PL.pl_id,1)
                                CSG.PL=CSG.PL.allocate(1);
                            end
                            if isfield(SH(ii,1),'DistMin');thresh.DistMin=SH(ii,1).DistMin;else thresh.DistMin=nan;end
                            if isfield(SH(ii,1),'DistMax');thresh.DistMax=SH(ii,1).DistMax;else thresh.DistMax=nan;end
                            if isfield(SH(ii,1),'LcMin');thresh.LcMin=SH(ii,1).LcMin;else thresh.LcMin=nan;end
                            if isfield(SH(ii,1),'LcMax');thresh.LcMax=SH(ii,1).LcMax;else thresh.LcMax=nan;end
                            CSG=CSG.addnewpoly(Xs{jj,1}',Ys{jj,1}',thresh);
                        end
                    end
                case 'Line'
                    CSG.dummy=2;
                    for ii=1:size(SH,1)
                        ii;
                        [Xs, Ys]=msim_polysplit(SH(ii,1).X,SH(ii,1).Y);
                        if isfield(SH(ii,1),'DistMin');thresh.DistMin=SH(ii,1).DistMin;else thresh.DistMin=nan;end
                        if isfield(SH(ii,1),'DistMax');thresh.DistMax=SH(ii,1).DistMax;else thresh.DistMax=nan;end
                        if isfield(SH(ii,1),'LcMin');thresh.LcMin=SH(ii,1).LcMin;else thresh.LcMin=nan;end
                        if isfield(SH(ii,1),'LcMax');thresh.LcMax=SH(ii,1).LcMax;else thresh.LcMax=nan;end
                        for jj=1:size(Xs,1)
                            for kk=1:size(Xs{jj,1},2)-1
                                x1=Xs{jj,1}(1,kk);y1=Ys{jj,1}(1,kk);
                                x2=Xs{jj,1}(1,kk+1);y2=Ys{jj,1}(1,kk+1);
                                CSG=CSG.addnewline([x1 y1 x2 y2],thresh);
                            end
                        end
                    end
                case 'Point'
                    CSG.dummy=3;
                    for ii=1:size(SH,1)
                        ii;
                        if isfield(SH(ii,1),'DistMin');thresh.DistMin=SH(ii,1).DistMin;else thresh.DistMin=nan;end
                        if isfield(SH(ii,1),'DistMax');thresh.DistMax=SH(ii,1).DistMax;else thresh.DistMax=nan;end
                        if isfield(SH(ii,1),'LcMin');thresh.LcMin=SH(ii,1).LcMin;else thresh.LcMin=nan;end
                        if isfield(SH(ii,1),'LcMax');thresh.LcMax=SH(ii,1).LcMax;else thresh.LcMax=nan;end
                        CSG=CSG.addpoint(SH(ii,1).X,SH(ii,1).Y,thresh);
                    end
            end
        end
        
        function CSG=addpoint(CSG,px,py,thresh)
            L_ID=[];
            %[L_ID P_ID]=CSG.GRD.getIDofNearFeatures_p(px,py);
            if CSG.LS.Nln>0%~isempty(L_ID)
                L_ID=1:CSG.LS.Nln;
                %check if the point intersects any line
                t=CSG.LS.isPointonLine(px,py,CSG.PS.coord(CSG.LS.l_id(L_ID,1),1),...
                                      CSG.PS.coord(CSG.LS.l_id(L_ID,1),2),...
                                      CSG.PS.coord(CSG.LS.l_id(L_ID,2),1),...
                                      CSG.PS.coord(CSG.LS.l_id(L_ID,2),2));
                id_cross_line=find(~isnan(t), 1);
                if ~isempty(id_cross_line)
                    p_cros=[CSG.PS.coord(CSG.LS.l_id(id_cross_line,1),1)+(CSG.PS.coord(CSG.LS.l_id(id_cross_line,2),1)-CSG.PS.coord(CSG.LS.l_id(id_cross_line,1),1))*t(id_cross_line) ...
                            CSG.PS.coord(CSG.LS.l_id(id_cross_line,1),2)+(CSG.PS.coord(CSG.LS.l_id(id_cross_line,2),2)-CSG.PS.coord(CSG.LS.l_id(id_cross_line,1),2))*t(id_cross_line)];
                    px=p_cros(1);py=p_cros(2);
                end
            end
            id=CSG.PS.isPointonSet([px py],[]);%P_ID
            if isnan(id)
                CSG.PS=CSG.PS.addnewpoint([px py],thresh);
                CSG.PS.pe=CSG.PS.Npnt;
                id=CSG.PS.pe;
                %CSG.GRD=CSG.GRD.assingPoint(px,py,id);
            else
                CSG.PS.pe=id;
                CSG.PS.DistMin(id)=nanmin(CSG.PS.DistMin(id),thresh.DistMin);
                CSG.PS.DistMax(id)=nanmax(CSG.PS.DistMax(id),thresh.DistMax);
                CSG.PS.LcMin(id)=nanmin(CSG.PS.LcMin(id),thresh.LcMin);
                CSG.PS.LcMax(id)=nanmax(CSG.PS.LcMax(id),thresh.LcMax);
            end
            CSG.PS.type(id)=CSG.dummy;
            if ~isempty(L_ID)
                if ~isempty(id_cross_line)
                    id_line1=CSG.LS.isLineinSet(CSG.LS.l_id(id_cross_line,1),id);
                    id_line2=CSG.LS.isLineinSet(id,CSG.LS.l_id(id_cross_line,2));
                    if ~((isinf(id_line1) && ~isnan(id_line2)) || (isinf(id_line2) && ~isnan(id_line1)))
                        %split the line
                        id_end=CSG.LS.l_id(id_cross_line,2);
                        CSG.LS.l_id(id_cross_line,2)=id;
                        temp=struct('DistMin',CSG.LS.DistMin(id_cross_line,1),...
                                    'DistMax',CSG.LS.DistMax(id_cross_line,1),...
                                    'LcMin',CSG.LS.LcMin(id_cross_line,1),...
                                    'LcMax',CSG.LS.LcMax(id_cross_line,1));
                        CSG.LS=CSG.LS.addline(id,id_end,temp);
                        id_new_line=CSG.LS.Nln;
                        %now we have split the line with id id_cross_line
                        %to two lines with id id_cross_line and id_new_line
                        %Now we have to find any feature that is made of
                        %the id_cross_line and update the structure
                        for ii=1:CSG.PL.Npl
                            for jj=1:length(CSG.PL.pl_id{ii})
                                if CSG.PL.pl_id{ii}(jj)==id_cross_line
                                    CSG.PL.pl_id{ii}=[CSG.PL.pl_id{ii}(1:jj);id_new_line;CSG.PL.pl_id{ii}(jj+1:end)];
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
        function CSG=addnewline(CSG,ln,thresh)
            %[L_ID P_ID]=CSG.GRD.getIDofNearFeatures_l(ln);
            %first check if the line intersects with other lines
            LN=ln;
            if CSG.LS.Nln>0%~isempty(L_ID)
                L_ID=1:CSG.LS.Nln;
                pint=CSG.LS.linexline(CSG.PS.coord(CSG.LS.l_id(L_ID,1),1), CSG.PS.coord(CSG.LS.l_id(L_ID,1),2), ...
                                      CSG.PS.coord(CSG.LS.l_id(L_ID,2),1), CSG.PS.coord(CSG.LS.l_id(L_ID,2),2), ...
                                  repmat(ln(1,1),length(L_ID),1),repmat(ln(1,2),length(L_ID),1),repmat(ln(1,3),length(L_ID),1),repmat(ln(1,4),length(L_ID),1));
                id_x=find(~isnan(pint(:,1)));
                if ~isempty(id_x)
                    t=CSG.LS.find_parametric_point(ln,pint(id_x,:));
                    [c_t, d_t]=sort(t,'ascend');
                    LN=[ln(1:2) pint(id_x(d_t(1)),:)];
                    for ii=2:length(c_t)
                        LN=[LN;pint(id_x(d_t(ii-1)),:) pint(id_x(d_t(ii)),:)];
                    end
                    LN=[LN;pint(id_x(d_t(end)),:) ln(3:4)];
                end
            end
            for ii=1:size(LN,1)
                CSG=CSG.addpoint(LN(ii,1),LN(ii,2),thresh);
                id1=CSG.PS.pe;
                CSG=CSG.addpoint(LN(ii,3),LN(ii,4),thresh);
                id2=CSG.PS.pe;
                %add line id1 id2
                %check if the line exist in the set
                id_line=CSG.LS.isLineinSet(id1, id2);
                if isinf(id_line)
                    continue
                end
                if isnan(id_line)
                    CSG.LS=CSG.LS.addline(id1,id2,thresh);
                    CSG.LS.le=CSG.LS.Nln;
                    %CSG.GRD=CSG.GRD.assingLine(ln,CSG.LS.le);
                else
                    CSG.LS.le=id_line;
                    CSG.LS.DistMin(abs(id_line))=nanmin(CSG.LS.DistMin(abs(id_line)),thresh.DistMin);
                    CSG.LS.DistMax(abs(id_line))=nanmax(CSG.LS.DistMax(abs(id_line)),thresh.DistMax);
                    CSG.LS.LcMin(abs(id_line))=nanmin(CSG.LS.LcMin(abs(id_line)),thresh.LcMin);
                    CSG.LS.LcMax(abs(id_line))=nanmax(CSG.LS.LcMax(abs(id_line)),thresh.LcMax);
                end
            end
        end
        
        function CSG=addnewpoly(CSG,xp,yp,thresh)
            CSG.PL.Npl=CSG.PL.Npl+1;
            for i=1:length(xp)-1
                CSG=addnewline(CSG,[xp(i) yp(i) xp(i+1) yp(i+1)],thresh);
                CSG.PL.pl_id{CSG.PL.Npl,1}=[CSG.PL.pl_id{CSG.PL.Npl,1};CSG.LS.le];
            end
        end
        
        function plotCSGobj(CSG,ln_cl,pnt_cl,lnw)
            if nargin<2
                ln_cl='b';
                pnt_cl='r';
                lnw=1.5;
            end
            hold on
            for i=1:CSG.LS.Nln
                plot(CSG.PS.coord(CSG.LS.l_id(i,:),1),CSG.PS.coord(CSG.LS.l_id(i,:),2),'color',ln_cl,'linewidth',lnw)
                drawnow
            end
            plot(CSG.PS.coord(:,1),CSG.PS.coord(:,2),'.','color',pnt_cl);
            
        end
        
        function writegeo(CSG,filename,opt)
            % writegeo(filename,meshoptions)
            fid=fopen([filename '.geo'],'w');
            fprintf(fid,'// ---- Mesh Size options ----//\n');
            fprintf(fid,'lc_gen = %f;\n',opt.lc_gen);
            
            %Group thresholds
            fields=[(1:CSG.PS.Npnt)' CSG.PS.DistMin(1:CSG.PS.Npnt) CSG.PS.DistMax(1:CSG.PS.Npnt)  CSG.PS.LcMin(1:CSG.PS.Npnt)  CSG.PS.LcMax(1:CSG.PS.Npnt) ; ...
                    (1:CSG.LS.Nln)' CSG.LS.DistMin(1:CSG.LS.Nln) CSG.LS.DistMax(1:CSG.LS.Nln) CSG.LS.LcMin(1:CSG.LS.Nln) CSG.LS.LcMax(1:CSG.LS.Nln)];
            id_fields=(1:size(fields,1))';
            id_nan=isnan(fields(:,2));
            id_fields(id_nan,:)=[];%fields(id_nan,:)=[] ;  
            [b, m, n]=unique(fields(id_fields,2:5),'rows');
            id_atr_thr=[(1:size(b,1))' size(b,1)+(1:size(b,1))'];
            for i=1:size(b,1)
                fprintf(fid,'b%g_DistMin = %f;\n',i,b(i,1));
                fprintf(fid,'b%g_DistMax = %f;\n',i,b(i,2));
                fprintf(fid,'b%g_LcMin = %f;\n',i,b(i,3));
                fprintf(fid,'b%g_LcMax = %f;\n',i,b(i,4));
            end
            fprintf(fid,'// ---Point Coordinates--------//\n');
            % find which points have minimum element size definition
            id_min_elem_def=find(~isnan(CSG.PS.LcMin(1:CSG.PS.Npnt)));
            id_no_def=find(isnan(CSG.PS.LcMin(1:CSG.PS.Npnt)));
            
            
            %write the points without definitions
            fprintf(fid,'Point(%i) = {%f, %f, %f, lc_gen};\n',[id_no_def  ...
                CSG.PS.coord(id_no_def,1) CSG.PS.coord(id_no_def,2)...
                zeros(length(id_no_def),1)]');
            %write the points with definitions
            for k=1:size(b,1)
                %find the id of the points associated with the i definition
                id_b= CSG.PS.LcMin(id_min_elem_def)==b(k,1);
                fprintf(fid,'Point(%i) = {%f, %f, %f, b%g_LcMin};\n',[id_min_elem_def(id_b)  ...
                        CSG.PS.coord(id_min_elem_def(id_b),1) CSG.PS.coord(id_min_elem_def(id_b),2)...
                        zeros(length(id_min_elem_def(id_b)),1) k*ones(length(id_min_elem_def(id_b)),1)]');
            end
            
            fprintf(fid,'// ---Lines indices--------//\n');
            %write the lines
            fprintf(fid,'Line(%i) = {%i,%i};\n',[(1:CSG.LS.Nln)',...
                CSG.LS.l_id(1:CSG.LS.Nln,1),CSG.LS.l_id(1:CSG.LS.Nln,2)]');
            
            %write line loops
            %It is always assumed that the first line loop corresponds to
            %domain outline
            fprintf(fid,'// ---Polygons--------//\n');
            for kk = 1:CSG.PL.Npl
                fprintf(fid,'Line Loop(%i) = {',kk);
                for ii = 1:length(CSG.PL.pl_id{kk})-1
                    fprintf(fid,'%i,',CSG.PL.pl_id{kk}(ii));
                end
                fprintf(fid,'%i};\n',CSG.PL.pl_id{kk}(ii+1));
            end
            
            %Write the Plane Surface
            fprintf(fid,'Plane Surface(%i) = {',1);
            for kk=1:CSG.PL.Npl-1
                fprintf(fid,'%i,',kk);
            end
            fprintf(fid,'%i};\n',CSG.PL.Npl);
            
            if opt.embed_points
                %embed points into Plane
                fprintf(fid,'Point{%i} In Surface{%i};\n',[(1:CSG.PS.Npnt)' ones(CSG.PS.Npnt,1)]');
            end
            if opt.embed_lines
                %embed lines into Plane
                fprintf(fid,'Line{%i} In Surface{%i};\n',[(1:CSG.LS.Nln)' ones(CSG.LS.Nln,1)]');
            end
            fprintf(fid,'// ---Attractors and Thresholds--------//\n');
            %Write fields
            %for each field write a list of nodes and edges
            for k=1:size(b,1)
                id=find(n==k);%find the thresholds associated with the b==k
                %find which are nodes and which are edges
                id_nd=find(id_fields(id)<=CSG.PS.Npnt);
                id_dg=find(id_fields(id)>CSG.PS.Npnt);
                
                fprintf(fid,'Field[%g] = Attractor;\n',id_atr_thr(k,1));
                if ~isempty(id_nd)
                    fprintf(fid,'Field[%g].NodesList = {',id_atr_thr(k,1));
                    for kk=1:length(id_nd)-1
                        fprintf(fid,'%i,',fields(id_fields(id(id_nd(kk))),1));
                    end
                    fprintf(fid,'%g};\n',fields(id_fields(id(id_nd(end))),1));
                end
                if ~isempty(id_dg)
                    fprintf(fid,'Field[%g].NNodesByEdge = %f;\n',id_atr_thr(k,1),b(k,1));
                    fprintf(fid,'Field[%g].EdgesList = {',id_atr_thr(k,1));
                    for kk=1:length(id_dg)-1
                        fprintf(fid,'%i,',fields(id_fields(id(id_dg(kk))),1));
                    end
                    fprintf(fid,'%i};\n',fields(id_fields(id(id_dg(end))),1));
                end
                
                fprintf(fid,'Field[%g] = Threshold;\n',id_atr_thr(k,2));
                fprintf(fid,'Field[%g].IField = %i;\n',id_atr_thr(k,2),id_atr_thr(k,1));
                fprintf(fid,'Field[%g].LcMin = b%g_LcMin;\n',id_atr_thr(k,2),k);
                fprintf(fid,'Field[%g].LcMax = b%g_LcMax;\n',id_atr_thr(k,2),k);
                fprintf(fid,'Field[%g].DistMin = b%g_DistMin;\n',id_atr_thr(k,2),k);
                fprintf(fid,'Field[%g].DistMax = b%g_DistMax;\n',id_atr_thr(k,2),k);
            end
            if size(b,1)>=1
                fprintf(fid,'Field[%g] = Min;\n',id_atr_thr(k,end)+1);
                fprintf(fid,'Field[%g].FieldsList = {',id_atr_thr(k,end)+1);
                for j=1:size(b,1)-1
                    fprintf(fid,'%i,',id_atr_thr(j,2));
                end
                fprintf(fid,'%i};\n',id_atr_thr(end,2));

                fprintf(fid,'Background Field = %i;\n',id_atr_thr(k,end)+1);
            end
            
            fprintf(fid,'// ---- Other Mesh options ----//\n');
            
            if strcmp(opt.el_type,'quad')
                fprintf(fid,'Mesh.SubdivisionAlgorithm = 1;\n');
                fprintf(fid,'Mesh.RecombineAll = 1;\n');
                fprintf(fid,'Mesh.Algorithm = 8;\n');
            end
            fprintf(fid,'Mesh.ElementOrder = %d;\n',opt.el_order);
            fprintf(fid,'Mesh.CharacteristicLengthExtendFromBoundary = 0;\n');
            fprintf(fid,'Mesh.CharacteristicLengthFromPoints = 0;\n');
            fprintf(fid,'Mesh.CharacteristicLengthFromCurvature = 0;\n');
            fprintf(fid,'Mesh.CharacteristicLengthMax = %f;\n',opt.lc_gen);
            fprintf(fid,'Mesh.SecondOrderIncomplete = %d;\n',opt.incomplete);
            
            %fprintf(fid,'// ---- Partition options ----//\n');
            %part_method=1;%1:chaco, 2:Metis
            %fprintf(fid,'Mesh.Partitioner = %g;\n',part_method);
            %fprintf(fid,'//Mesh.NbPartitions = %g;\n',7);
            %fprintf(fid,'Mesh.ColorCarousel = %g;\n',3);
            %fprintf(fid,'Mesh.MshFilePartitioned = 0;\n');
            %if part_method==1
            %    fprintf(fid,'Mesh.ChacoGlobalMethod = 2;\n');
            %    fprintf(fid,'Mesh.ChacoArchitecture = 1;\n');
            %    fprintf(fid,'Mesh.ChacoParamREFINE_MAP = 0;\n');
            %else
            %    fprintf(fid,'Mesh.MetisAlgorithm = 1;\n');
            %end
            fprintf(fid,'// ------------------------//\n');
            
            fclose(fid);
        end
        
        function runGmsh(CSG,filename,gmsh_exe, Nsub)
            if isempty(gmsh_exe)
                gmsh_exe='E:\downloads\gmsh-2.6.1-Windows\gmsh.exe';
            end
            if isempty(Nsub)
                system([gmsh_exe ' ' filename '.geo -2']);
            else
                system([gmsh_exe ' ' filename '.geo -2 -part ' num2str(Nsub)]);
            end
        end
        
        function showGmsh(CSG,filename,gmsh_exe,type)
            if isempty(gmsh_exe)
                gmsh_exe='E:\downloads\gmsh-2.6.1-Windows\gmsh.exe';
            end
            switch type
                case 'geo'
                    system([gmsh_exe ' ' filename '.geo']);
                case 'msh'
                    system([gmsh_exe ' ' filename '.msh']);
                case  'geomsh'
                    system([gmsh_exe ' ' filename '.geo ' filename '.msh']);
            end
        end
        
        function CSG=readShapefiles_withArc(CSG,nonpointlist,pointlist,tempname)
            warning('This function has not been thoroughly tested')
            
            button = questdlg({['WARNING: This script will DELETE ALL FILES STARTING WITH'],...
                                [' '],...
                                [tempname '*']});
                    if strcmp(button,'Cancel') || strcmp(button,'No')
                        return
                    end
            
            python_cmnd='''createCSG_v2.py''';
            python_cmnd=[python_cmnd ',' '''' nonpointlist{1,1} ''''];
            for i=2:size(nonpointlist,1)
                python_cmnd=[python_cmnd ',' '''' nonpointlist{i,1} ''''];
            end
            python_cmnd=[python_cmnd ',' '''' tempname '''' ',' '''' num2str(CSG.tol) ''''];
            %delete all files with names starting from "tempname"
            delete([tempname '*'])
            eval(['msim_python(' python_cmnd ')']);
            
            %make a list of points in CSG.
            Verts=shaperead([tempname '_allVert']);
            for i=1:size(Verts,1);
                id=CSG.PS.isPointonSet([Verts(i,1).X Verts(i,1).Y],[]);
                if ~isnan(id);continue;end
                i;
                if ~isfield(Verts(i,1),'LcMin')
                    thresh.DistMin=nan; thresh.DistMax=nan;
                    thresh.LcMin=nan; thresh.LcMax=nan;
                elseif Verts(i,1).LcMin==0;
                    thresh.DistMin=nan; thresh.DistMax=nan;
                    thresh.LcMin=nan; thresh.LcMax=nan;
                else
                    thresh.DistMin=Verts(i,1).DistMin; thresh.DistMax=Verts(i,1).DistMax;
                    thresh.LcMin=Verts(i,1).LcMin; thresh.LcMax=Verts(i,1).LcMax;
                end
                CSG.PS=CSG.PS.addnewpoint([Verts(i,1).X Verts(i,1).Y],thresh);
            end
            
            %find out which vertices compose the somain outline
            %we always assume that only one polygon is the domain and the others are
            %just holes. The domain has clockwise orientation
            dom=shaperead([tempname '_poly_diss_int']);
            id_nan=find(isnan(dom(1,1).X));
            N_poly=length(id_nan);
            id_nan=[0 id_nan];
            poly_orient=nan(N_poly,1);
            for i=1:N_poly
                poly_orient(i)=ispolycw(dom(1,1).X(id_nan(i)+1:id_nan(i+1)-1), dom(1,1).Y(id_nan(i)+1:id_nan(i+1)-1))==1;
            end
            if sum(poly_orient)>1
                error('More than one polygons have counter clock wise orientation')
            end
            [cc dd]=sort(poly_orient,'descend');
            
            cnt=1;
            for i=1:N_poly
                poly_id_start(i)=cnt;
                x_p=dom(1,1).X(id_nan(dd(i))+1:id_nan(dd(i)+1)-1);
                y_p=dom(1,1).Y(id_nan(dd(i))+1:id_nan(dd(i)+1)-1);
                for j=1:length(x_p)-1;
                    id1=CSG.PS.isPointonSet([x_p(j) y_p(j)],[]);
                    if isnan(id1)
                        if ~isfield(dom(1,1),'LcMin')
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                        elseif Verts(i,1).LcMin==0;
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                        else
                            thresh.DistMin=dom(1,1).DistMin; thresh.DistMax=dom(1,1).DistMax;
                            thresh.LcMin=dom(1,1).LcMin; thresh.LcMax=dom(1,1).LcMax;
                        end
                        CSG.PS=CSG.PS.addnewpoint([x_p(j) y_p(j)],thresh);
                        id1=CSG.PS.Npnt;
                    end
                    id2=CSG.PS.isPointonSet([x_p(j+1) y_p(j+1)],[]);
                    if isnan(id2)
                        if ~isfield(dom(1,1),'LcMin')
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                        elseif dom(1,1).LcMin==0;
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                        else
                            thresh.DistMin=dom(1,1).DistMin; thresh.DistMax=dom(1,1).DistMax;
                            thresh.LcMin=dom(1,1).LcMin; thresh.LcMax=dom(1,1).LcMax;
                        end
                        CSG.PS=CSG.PS.addnewpoint([x_p(j+1) y_p(j+1)],thresh);
                        id2=CSG.PS.Npnt;
                    end
                    idline=CSG.LS.isLineinSet(id1, id2);
                    if isnan(idline)
                        thresh.DistMin=nanmin(CSG.PS.DistMin(id1),CSG.PS.DistMin(id2));
                        thresh.DistMax=nanmin(CSG.PS.DistMax(id1),CSG.PS.DistMax(id2));
                        thresh.LcMin=nanmin(CSG.PS.LcMin(id1),CSG.PS.LcMin(id2));
                        thresh.LcMax=nanmin(CSG.PS.LcMax(id1),CSG.PS.LcMax(id2));
                        CSG.LS=CSG.LS.addline(id1,id2,thresh);
                        cnt=cnt+1;
                    end
                end
                poly_id_end(i)=cnt-1;
                CSG.PL.pl_id{i,1}=poly_id_start(i):poly_id_end(i);
                CSG.PL.Npl=CSG.PL.Npl+1;
            end
            
            %add the remaining lines
            lines=shaperead(tempname);
            for k=1:size(lines,1)
                id_nan=find(isnan(lines(k,1).X));
                N_lines=length(id_nan);
                id_nan=[0 id_nan];
                for i=1:N_lines
                    [k i];
                    x_l=lines(k,1).X(id_nan(i)+1:id_nan(i+1)-1);
                    y_l=lines(k,1).Y(id_nan(i)+1:id_nan(i+1)-1);
                    for j=1:length(x_l)-1;
                        id1=CSG.PS.isPointonSet([x_l(j) y_l(j)],[]);
                        if isnan(id1)
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                            CSG.PS=CSG.PS.addnewpoint([x_l(j) x_l(j)],thresh);
                            id1=CSG.PS.Npnt;
                        end
                        id2=CSG.PS.isPointonSet([x_l(j+1) y_l(j+1)],[]);
                        if isnan(id2)
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                            CSG.PS=CSG.PS.addnewpoint([x_l(j+1) x_l(j+1)],thresh);
                            id2=CSG.PS.Npnt;
                        end
                        idline=CSG.LS.isLineinSet(id1, id2);
                        if isnan(idline)
                            thresh.DistMin=nanmin(CSG.PS.DistMin(id1),CSG.PS.DistMin(id2));
                            thresh.DistMax=nanmin(CSG.PS.DistMax(id1),CSG.PS.DistMax(id2));
                            thresh.LcMin=nanmin(CSG.PS.LcMin(id1),CSG.PS.LcMin(id2));
                            thresh.LcMax=nanmin(CSG.PS.LcMax(id1),CSG.PS.LcMax(id2));
                            CSG.LS=CSG.LS.addline(id1,id2,thresh);
                        end
                    end
                end
            end
            
            %last add the point features
            for i=1:size(pointlist,1)
                pnt=shaperead(pointlist{i,1});
                for j=1:size(pnt,1)
                    [i j];
                    id=CSG.PS.isPointonSet([pnt(j,1).X pnt(j,1).Y],[]);
                    if isnan(id)
                        if ~isfield(pnt(j,1),'LcMin')
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                        elseif pnt(j,1).LcMin==0;
                            thresh.DistMin=nan; thresh.DistMax=nan;
                            thresh.LcMin=nan; thresh.LcMax=nan;
                        else
                            thresh.DistMin=pnt(j,1).DistMin; thresh.DistMax=pnt(j,1).DistMax;
                            thresh.LcMin=pnt(j,1).LcMin; thresh.LcMax=pnt(j,1).LcMax;
                        end
                        CSG.PS=CSG.PS.addnewpoint([pnt(j,1).X pnt(j,1).Y],thresh);
                    end
                end
            end
            
        end
    end
end
