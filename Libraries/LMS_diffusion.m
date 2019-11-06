function [ wyj ] = LMS_diffusion( alpha, a1, specificHeat, k, rho, Tmax, dt, U1)
    U=U1;
    a0=sqrt(rho*specificHeat/k);
    UT=Fractional_diffusion_solution( alpha, 500/(a0)^2, a0, Tmax, dt)*a1;
    ER=0;
    for kk=1:Tmax/dt
        ER=ER+(U(kk)-UT(kk))^2;
    end
wyj=ER;
end
