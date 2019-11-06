function [u,pb]=FEM_plate_sim(c,a,f,d,geometry,initial_u,t,bc_fun,mash_max_size,exact_geometry_match,pb0,linear)

    numberOfPDE = 1;
    pb = createpde(numberOfPDE);
    geometryFromEdges(pb,geometry);
    applyBoundaryCondition(pb,'neumann','Edge',2,'g',bc_fun,'q',0,'Vectorized','off');
    msh=generateMesh(pb,'Hmax',mash_max_size,'MesherVersion', 'R2013a','GeometricOrder','linear');
    
    if ~exact_geometry_match
        if linear
            R=RegBoundCond(initial_u(:,end), pb0.Mesh.Nodes);
            [xData, yData] = prepareCurveData( R(:,2), R(:,1) );

            ft = fittype( 'exp2' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

            [fitresult, gof] = fit( xData, yData, ft, opts );
            u0=fitresult(msh.Nodes(2,:));
        else
            results = createPDEResults(pb0,initial_u(:,end));
            u0=interpolateSolution(results,msh.Nodes(1,:),msh.Nodes(2,:)); 
        end
    else
        if size(initial_u)==size(msh.Nodes)
            u0=initial_u;
        else
            u0=zeros(size(msh.Nodes,2),1);
        end
    end
    
    u = parabolic(u0,t,pb,c,a,f,d);
end

% 
% 
% k=800;
% rho = 8960;
% specificHeat = 386;
% thick = .01;
% stefanBoltz = 5.670373e-8;
% hCoeff = 0;
% ta = 0;
% emiss = 0;
% width=1.5;
% height=1.5;
% hhh=1.5;
% ww=width/9;
% hh=hhh/9;
% gam=(delta0*width/6)^2;
% gam2=(deltaK*width/6)^2;
% c = sprintf('%g*%g', thick, k);%9999175
% c2 = sprintf('(-(floor(floor(1+%g-(x-0.75).^2-(y-0.75).^2)))*%g+1)*%g*%g', gam2, aM, thick, k);%9999175
% a = sprintf('2*%g + 2*%g*%g*u.^3', hCoeff, emiss, stefanBoltz);
% f = 2*hCoeff*ta + 2*emiss*stefanBoltz*ta^4;
% d = thick*rho*specificHeat;
% % c = 0.1;
% % a = 0;
% % f = 0;
% % d = 10;
% % Plot the geometry and display the edge labels for use in the boundary
% % condition definition.
% gdmTrans = [3 4 0  width width 0 0 0 height height];
% % Create a pde entity for a PDE with a single dependent variable
% numberOfPDE = 1;
% pb = createpde(numberOfPDE);
% pb0=createpde(numberOfPDE);
% % Create a geometry entity
% sf2=['R1'];
% ns2=[82; 49];
% g0 = geomFractal0;
% g = geomFractal;
% hmax0=hmax;
% %figure;
% %pg = pdeGeometryFromEdges(g);
% %pdegplot(g0, 'edgeLabels', 'on');
% % Solution is zero at all four outer edges of the square
% %pb.BoundaryConditions = pdeBoundaryConditions(pg.Edges(1), 'u', 0);
% %pb.BoundaryConditions = pdeBoundaryConditions(pg.Edges(1), 'u', 0.5);
% %fun=@(region,state)3*sin(state.time*0.5);
% fun=@(region,state)500*k*thick;%*sin(state.time);
% BC = applyBoundaryCondition(pb,'neumann','Edge',2,'g',fun,'q',0,'Vectorized','off');
% applyBoundaryCondition(pb,'neumann','Edge',2,'g',fun,'q',0,'Vectorized','off');
% % [p, e, t] = initmesh(g0, 'Hmax', hmax0);
% %[p2, e2, t2] = initmesh(g, 'Hmax', hmax);
% geometryFromEdges(pb,g)
% geometryFromEdges(pb0,g0)
% applyBoundaryCondition(pb0,'neumann','Edge',2,'g',fun,'q',0,'Vectorized','off');
% 
% msh=generateMesh(pb,'Hmax',hmax,'MesherVersion', 'R2013a','GeometricOrder','linear');
% msh0=generateMesh(pb0,'Hmax',hmax,'MesherVersion', 'R2013a','GeometricOrder','linear');
% 
% u0 = zeros(size(msh0.Nodes,2),1);
% %ix = find(sqrt(p(1,:).^2+p(2,:).^2)<0.4);
% %u0(ix) = ones(size(ix));
% 
% nframes = 1000;
% tlist = linspace(0,1000,nframes);
% 
% geometryFromEdges(pb0,g0)
% u1 = parabolic(u0,tlist,pb0,c,a,f,d);
% 
% 
% results = createPDEResults(pb0,u1,tlist,'time-dependent');
% 
% results = createPDEResults(pb0,u1(:,end));
% 
% u0=interpolateSolution(results,msh.Nodes(1,:),msh.Nodes(2,:)); 
% % R=RegBoundCond(u1(:,end), msh0.Nodes);
% % [xData, yData] = prepareCurveData( R(:,2), R(:,1) );
% % 
% % ft = fittype( 'exp2' );
% % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % 
% % [fitresult, gof] = fit( xData, yData, ft, opts );
% % 
% % u0=fitresult(msh.Nodes(2,:));
% 
% 
% nframes = 1000;
% 
% tlist2 = linspace(0,1000,nframes);
% 
% u3 = parabolic(u0,tlist,pb,c,a,f,d);
% 
% % figure;
% %   pdeplot(pb0, 'xydata', u1(:,end)/10, 'contour', 'on', 'colormap', 'jet');
% %    title(sprintf('Temperatura płytki, analiza czasowa (%d sekund)\n', ...
% %     tlist(1,end)));
% % 
% % figure;
% %   pdeplot(pb, 'xydata', u0(:,1)/10, 'contour', 'on', 'colormap', 'jet');
% %    title(sprintf('Temperatura płytki, analiza czasowa (%d sekund)\n', ...
% %     tlist(1,end)));
% % figure;
% %   pdeplot(pb, 'xydata', u3(:,end)/10, 'contour', 'on', 'colormap', 'jet');
% %    title(sprintf('Temperatura płytki, analiza czasowa (%d sekund)\n', ...
% %     tlist(1,end)));
% % for tt=1:size(p,2)
% %     if p(1,tt)<(0.75+hmax0) && p(1,tt)>(0.75-hmax0) && p(2,tt)==0
% %         
% %         u2=u1(tt,:);
% %     end
% % end
% % u2TMP=[];
% % for tt=1:size(msh.Nodes,2)
% %     if msh.Nodes(2,tt)==0
% %         u2TMP=[u2TMP;u3(tt,:)];
% %     end
% % end
% size(u0(:,1))
% size(u1(:,end))
% 
% E1=sum(u0(:,1))/area(pb.Mesh);
% E2=sum(u1(:,end))/area(pb0.Mesh);
% 
% E1=mean(u0(:,1))*area(pb.Mesh);
% E2=mean(u1(:,end))*area(pb0.Mesh);
% 
% end

