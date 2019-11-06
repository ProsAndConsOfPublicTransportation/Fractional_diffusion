function wyj=Fractional_diffusion_solution( alpha, U0, a, tmax, dt)
    for t=dt:dt:tmax
        UW(t)=U0*2*a/(alpha*gamma(alpha/2))*(t-dt)^(alpha/2);
    end
    wyj=UW;
end
