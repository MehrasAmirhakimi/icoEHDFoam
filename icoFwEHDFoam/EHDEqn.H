        solve
        (
            fvm::laplacian(phii) + chargeDensity/epsilon
        );
        posFlux = -bpos*mesh.magSf()*fvc::snGrad(phii);
        negFlux = -bneg*mesh.magSf()*fvc::snGrad(phii);
        
        elec=-fvc::grad(phii);
		omega=sqrt(pow(electron,3)/(4*pinum*epsilon*pow(kB*T0,2))*mag(elec)) ; 

        forAll(fw,celli)
		{
			twoOmega=2*omega[celli];

			
			if (abs(twoOmega)<scalar(3.75))
			{
				yy=pow(twoOmega/3.75,2);
				bessi1=twoOmega*(p1+yy*(p2+yy*(p3+yy*(p4+yy*(p5+yy*(p6+yy*p7))))));
			}
			else
			{
				yy=3.75/abs(twoOmega);
				ax=abs(twoOmega);
				bessi1=(Foam::exp(ax)/Foam::sqrt(ax))*(q1+yy*(q2+yy*(q3+yy*(q4+yy*(q5+yy*(q6+yy*(q7+yy*(q8+yy*q9))))))));
				if (twoOmega<0)
				{
					bessi1=-bessi1;
				}
			}
			if (twoOmega==0)
			{
				fw[celli]=1.0;
			}
			else 
			{
				fw[celli]=bessi1/(twoOmega/2);
			}

		}

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
