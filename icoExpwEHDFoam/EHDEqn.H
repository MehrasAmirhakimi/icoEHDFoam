        solve
        (
            fvm::laplacian(phii) + chargeDensity/epsilon
        );
        posFlux = -bpos*mesh.magSf()*fvc::snGrad(phii);
        negFlux = -bneg*mesh.magSf()*fvc::snGrad(phii);
        
        elec=-fvc::grad(phii);
		omega= 2*sqrt(pow(electron,3)*mag(elec)/epsilon)/(kB*T0);
		fw = Foam::exp(omega);

		source = fw*kdN - kr*pos*neg;
        solve
        (
            fvm::ddt(pos) 
            + fvm::div(posFlux, pos) 
            + fvm::div(phi, pos)             
            - fvm::laplacian(D,pos) 
            == source
        );
        solve
        (
            fvm::ddt(neg) 
            + fvm::div(negFlux, neg) 
            + fvm::div(phi, neg)              
            - fvm::laplacian(D,neg)
            == source
        );
        chargeDensity = pos - neg;
		force = -chargeDensity*fvc::grad(phii)/rho0;
