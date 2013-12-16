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
		// Reasonable parameters to nondimensionalize by.
		GLFloat N_QG = 1.3; // cm
		GLFloat T_QG = 12; // days
		GLFloat L_QG = 47; // km
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:256 domainMin:-1500/L_QG length:2000/L_QG];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:128 domainMin:-500/L_QG length:1000/L_QG];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		NSArray *spectralDimensions = [x dimensionsTransformedToBasis: x.differentiationBasis];
		
		GLLinearTransform *laplacian = [GLLinearTransform harmonicOperatorFromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *laplacianMinusOne = [laplacian plus: @(-1.0)];
		GLLinearTransform *inverseLaplacianMinusOne = [laplacianMinusOne inverse];
		
		GLLinearTransform *diff_xxx = [GLLinearTransform differentialOperatorWithDerivatives:@[@(3),@(0)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *diff_xyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(2)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *diff_xxy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(2),@(1)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *diff_yyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(0),@(3)] fromDimensions:spectralDimensions forEquation:equation];
		
		GLLinearTransform *diffJacobianX = [diff_xxx plus: diff_xyy];
		GLLinearTransform *diffJacobianY = [diff_xxy plus: diff_yyy];
		
		GLFloat k = 0.5*xDim.sampleInterval;
		GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *svv = [GLLinearTransform spectralVanishingViscosityFilterWithDimensions: spectralDimensions scaledForAntialiasing: YES forEquation: equation];
		GLLinearTransform *diffLin = [[biharmonic times: @(k)] times: svv];
		
		
		/************************************************************************************************/
		/*		Create the initial conditions															*/
		/************************************************************************************************/
		
		// We're going to plop down a gaussian = amplitude * exp( - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(length*length) );
		GLFloat amplitude = 15.0/N_QG;
		GLFloat length = 80/L_QG;
		
		GLFunction *r2 = [[x times: x] plus: [y times: y]];
		GLFunction *gaussian = [[[r2 times: @(-1.0/(length*length))] exponentiate] times: @(amplitude)];
		
		/************************************************************************************************/
		/*		Create a file to output data															*/
		/************************************************************************************************/
		
		// Now we create a mutable variable in order to record the evolution of the Gaussian.
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Desktop/Quasigeostrophy.nc"] forEquation: equation overwriteExisting: YES];
		GLMutableVariable *sshHistory = [gaussian variableByAddingDimension: tDim];
		sshHistory.name = @"SSH";
		sshHistory = [netcdfFile addVariable: sshHistory];
		
		/************************************************************************************************/
		/*		Estimate the time step size																*/
		/************************************************************************************************/
		
		GLFunction *v = [[gaussian x] spatialDomain];
		GLFunction *u = [[gaussian y] spatialDomain];
		GLFunction *speed = [[u times: u] plus: [v times: v]];
		[equation solveForVariable: speed];
		
		CGFloat cfl = 0.5;
		GLFloat U = sqrt([speed maxNow]);
		GLFloat timeStep = cfl * xDim.sampleInterval / U;
		
		/************************************************************************************************/
		/*		Create an integrator: dy/dt=f															*/
		/************************************************************************************************/
		
		GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[[gaussian differentiateWithOperator: laplacianMinusOne]] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
//		GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[[gaussian differentiateWithOperator: laplacianMinusOne]] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {

			
			// First, invert psi to get eta where (\nabla^2 -1) eta = psi. We use our cached differential operator.
			GLFunction *eta = [inverseLaplacianMinusOne transform: yNew[0]];
			
			// Second, compute f. For QG, f = (eta_{xxx} + eta_{xyy})*eta_y - (eta_{xxy}+eta_{yyy})*eta_x - eta_x + k*(eta_{xx}+eta_{yy})
			GLFunction *f = [[eta differentiateWithOperator: diffLin] plus: [[[[eta y] times: [eta differentiateWithOperator: diffJacobianX]] minus: [[eta x] times: [[[eta differentiateWithOperator: diffJacobianY] spatialDomain] plus: @(1.0)]]] frequencyDomain]];
			return @[f];
		}];
		
		/************************************************************************************************/
		/*		Now iterate! Stop every day to write out some data.										*/
		/************************************************************************************************/
		
		for (GLFloat time = 1/T_QG; time < 365/T_QG; time += 1/T_QG)
		{
            @autoreleasepool {
				y = [integrator stepForwardToTime: time][0];
				if (!y)	break;
				
				NSLog(@"Logging day: %f, last step size: %f, next step size: %f.", (time*T_QG), integrator.lastStepSize*T_QG, integrator.stepSize*T_QG);
				
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: @(time)];
				GLFunction *eta = [[inverseLaplacianMinusOne transform: y] spatialDomain];
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

