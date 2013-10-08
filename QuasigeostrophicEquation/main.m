//
//  main.m
//  Quasigeostrophy
//
//  Created by Jeffrey Early on 4/13/12.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main (int argc, const char * argv[])
{
	
	@autoreleasepool {
	    [GLVariable setPrefersSpatialMultiplication: YES];
        
		// Reasonable parameters to nondimensionalize by.
		GLFloat N_QG = 1.3; // cm
		GLFloat T_QG = 12; // days
		GLFloat L_QG = 47; // km
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		GLDimension *xDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 256 domainMin: -1500/L_QG length: 2000/L_QG];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 128 domainMin: -500/L_QG length: 1000/L_QG];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		// Now we choose the default differentiation basis (exponential, sine, or cosine);
		[equation setDefaultDifferentiationBasis: @[@(kGLSineHalfShiftBasis), @(kGLExponentialBasis)] forOrder: 2];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLVariable *y = [GLVariable variableOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		// Differential operators don't need to be cached like this, but in some cases it will give a performance boost to do so.
		// You could just write out the equation in the integrator.
		
		// At the moment we know that this is the spectral operators, although in the future we'll have to set this up explicitly.
		GLSpectralDifferentialOperatorPool *diffOperators = [equation defaultDifferentialOperatorPoolForVariable: x];
		
		// Create the operator xx+yy-1---this is how you compute y from eta
		GLSpectralDifferentialOperator *laplacianMinusOne = [[diffOperators harmonicOperator] scalarAdd: -1.0];
		[diffOperators setDifferentialOperator: laplacianMinusOne forName: @"laplacianMinusOne"];
		
		// Create the operator 1/(xx+yy-1)---this is how you compute eta from y.
		[diffOperators setDifferentialOperator: [laplacianMinusOne scalarDivide: 1.0] forName: @"inverseLaplacianMinusOne"];
		
		// This builds the differentiation matrix diff_{xxx} + diff_{xyy}
		[diffOperators setDifferentialOperator: [[diffOperators xxx] plus: [diffOperators xyy]] forName: @"diffJacobianX"];
		
		// This builds the differentiation matrix diff_{xxy} + diff_{yyy}
		[diffOperators setDifferentialOperator: [[diffOperators xxy] plus: [diffOperators yyy]] forName: @"diffJacobianY"];
		
		GLFloat k = 0.5*xDim.sampleInterval;
		GLSpectralDifferentialOperator *svv = [diffOperators spectralVanishingViscosityFilter];
		//		[diffOperators setDifferentialOperator: [[[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: k] minus: [diffOperators x]] forName: @"diffLin"];
		//		[diffOperators setDifferentialOperator: [[[[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: k] multiply: svv] minus: [diffOperators x]] forName: @"diffLin"];
		[diffOperators setDifferentialOperator: [[[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: k] multiply: svv] forName: @"diffLin"];
		
		
		/************************************************************************************************/
		/*		Create the initial conditions															*/
		/************************************************************************************************/
		
		// We're going to plop down a gaussian = amplitude * exp( - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(length*length) );
		GLFloat amplitude = 15.0/N_QG;
		GLFloat length = 80/L_QG;
		
		GLVariable *r2 = [[x times: x] plus: [y times: y]];
		GLVariable *gaussian = [[[r2 scalarMultiply: -1.0/(length*length)] exponentiate] scalarMultiply: amplitude];
		
		/************************************************************************************************/
		/*		Create a file to output data															*/
		/************************************************************************************************/
		
		// Now we create a mutable variable in order to record the evolution of the Gaussian.
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Desktop/Quasigeostrophy_pos.nc"] forEquation: equation overwriteExisting: YES];
		GLMutableVariable *sshHistory = [gaussian variableByAddingDimension: tDim];
		sshHistory.name = @"SSH";
		sshHistory = [netcdfFile addVariable: sshHistory];
		
		/************************************************************************************************/
		/*		Estimate the time step size																*/
		/************************************************************************************************/
		
		GLVariable *v = [[gaussian x] spatialDomain];
		GLVariable *u = [[gaussian y] spatialDomain];
		GLVariable *speed = [[u times: u] plus: [v times: v]];
		[equation solveForVariable: speed];
		
		CGFloat cfl = 0.5;
		GLFloat U = sqrt([speed maxNow]);
		GLFloat timeStep = cfl * xDim.sampleInterval / U;
		
		/************************************************************************************************/
		/*		Create an integrator: dy/dt=f															*/
		/************************************************************************************************/
		
		GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[[gaussian diff: @"laplacianMinusOne"]] stepSize: timeStep fFromTY:^(GLVariable *time, NSArray *yNew) {
			
			// First, invert psi to get eta where (\nabla^2 -1) eta = psi. We use our cached differential operator.
			GLVariable *eta = [yNew[0] diff: @"inverseLaplacianMinusOne"];
			
			// Second, compute f. For QG, f = (eta_{xxx} + eta_{xyy})*eta_y - (eta_{xxy}+eta_{yyy})*eta_x - eta_x + k*(eta_{xx}+eta_{yy})
			// Notice again that we're using the cached differential operators from above to do this.
			//GLVariable *f = [[eta diff:@"diffLin"] plus: [[[[eta y] times: [eta diff: @"diffJacobianX"]] minus: [[eta x] times: [eta diff: @"diffJacobianY"]]] frequencyDomain]];
			GLVariable *f = [[eta diff:@"diffLin"] plus: [[[[eta y] times: [eta diff: @"diffJacobianX"]] minus: [[eta x] times: [[[eta diff: @"diffJacobianY"] spatialDomain] scalarAdd: 1.0]]] frequencyDomain]];
			return @[f];
		}];
		
		/************************************************************************************************/
		/*		Now iterate! Stop every day to write out some data.										*/
		/************************************************************************************************/
		
		for (GLFloat time = 1/T_QG; time < 800/T_QG; time += 1/T_QG)
		{
            @autoreleasepool {
				y = [integrator stepForwardToTime: time][0];
				if (!y)	break;
				
				NSLog(@"Logging day: %f, last step size: %f, next step size: %f.", (integrator.currentTime*T_QG), integrator.lastStepSize*T_QG, integrator.stepSize*T_QG);
				
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: @(integrator.currentTime)];
				GLVariable *eta = [[y diff: @"inverseLaplacianMinusOne"] spatialDomain];
				[sshHistory concatenateWithLowerDimensionalVariable: eta alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
            }
		}
		
		NSLog(@"Close the NetCDF file and wrap up");
		
		[equation waitUntilAllOperationsAreFinished];
		
		// The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
		[netcdfFile waitUntilAllOperationsAreFinished];
		[netcdfFile close];
	}
    return 0;
}

