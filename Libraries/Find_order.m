function [ x,fval ] = Find_order(specificHeat,k,rho,u)
    funkcja=@(x) LMS_diffusion(x(1),x(2),specificHeat,k,rho,1000,1,u);
    [x,fval]=fminsearch(funkcja, [1, 1]);
end
