//Function to optimize MRF energy using TRW-S or BP algorithm (wrapper to Vladimir Kolmogorov's code).
//This version assumes that pairsise potentials can be "decomposed": V_{ij}(k,l) = P(i, j) * M(k, l).
//Here P depends only on variable indeces, M depends only on labels.
//Input examples:
//mrfMinimizeMex(U, P)
//mrfMinimizeMex(U, P, M)
//mrfMinimizeMex(U, P, M, options)
//Output examples:
//S = mrfMinimizeMex(U, P, M, options)
//[S, E] = mrfMinimizeMex(U, P, M, options)
//[S, E, LB] = mrfMinimizeMex(U, P, M, options)
//
//INPUT:
//	U		- unary terms (double[numLabels, numNodes])
//	P		- matrix of edge coefficients (sparse double[numNodes, numNodes]);
//	M		- matrix of label dependencies (double[naumLabels, numLabels]); if M is not specified, Potts is assumed
//				if you want to set options without M call: mrfMinimizeMex(U, P, [], options)
//options	- Stucture that determines metod to be used.
//				Fields:  
//					method		:	method to use (string: 'trw-s' or 'bp') default: 'trw-s'
//					maxIter		:	maximum number of iterations (double) default: 100
//					funcEps		:	If functional change is less than funcEps then stop, TRW-S only (double) default: 1e-2
//					verbosity	:	verbosity level: 0 - no output; 1 - final output; 2 - full output (double) default: 0
//					printMinIter:	After printMinIter iterations start printing the lower bound (double) default: 10
//					printIter	:	and print every printIter iterations (double) default: 5
//
//OUTPUT: 
//	S		- Labeling that has energy E, vector numNodes * 1 of type double (indeces are in [1,...,numLabels])
// 	E		- Energy #iter * 1 matrix of type double
// 	LB		- Lower bound #iter * 1 matrix of type double (only for TRW-S method)
// 
//  by Anton Osokin (firstname.lastname@gmail.com), 2011


#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <time.h>

#include "MRFEnergy.h"
#include "mex.h"

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { mexErrMsgTxt(msg);}

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
typedef int mwSize;
typedef int mwIndex;
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MATLAB_ASSERT(nrhs >= 2 , "Not enough input arguments, expected 2 - 4" ); \
	MATLAB_ASSERT(nrhs <= 4, "Too many input arguments, expected 2 - 4"); 	
	
	//Fix input parameter order:
	const mxArray *uInPtr = (nrhs > 0) ? prhs[0] : NULL; //unary
	const mxArray *pInPtr = (nrhs > 1) ? prhs[1] : NULL; //pairwise
	const mxArray *mInPtr = (nrhs > 2) ? prhs[2] : NULL; //label matrix
	const mxArray *oInPtr = (nrhs > 3) ? prhs[3] : NULL; //options
	
	//Fix output parameter order:
	mxArray **eOutPtr = (nlhs > 1) ? &plhs[1] : NULL; //energy
	mxArray **sOutPtr = (nlhs > 0) ? &plhs[0] : NULL; //solution
	mxArray **lbOutPtr = (nlhs > 2) ? &plhs[2] : NULL; //lowerbound
	mxArray **tOutPtr = (nlhs > 3) ? &plhs[3] : NULL; //time array

	//prepare default options
	MRFEnergy<TypeGeneral>::Options options;
	options.m_eps = 1e-2; 
	options.m_iterMax = 20;
	options.m_printIter = 5;     
	options.m_printMinIter = 10;
	int verbosityLevel = 0;
	int method = 0;
	
	//get options structure
	if(oInPtr != NULL){
		MATLAB_ASSERT(mxIsStruct(oInPtr), "Expected structure array for options");
		MATLAB_ASSERT(mxGetNumberOfElements(oInPtr) == 1, "Wrong size of options structure: expected 1");
		mxArray *curField = NULL;
		if((curField = mxGetField(oInPtr, 0, "method")) != NULL){
			MATLAB_ASSERT(mxGetClassID(curField) == mxCHAR_CLASS, "Wrong structure type for options: expected STRING for field <<method>>");
		
			mwSize buflen = mxGetN(curField)*sizeof(mxChar)+1;
			char *buf = (char*)mxMalloc(buflen);
			if(!mxGetString(curField, buf, buflen)){
				if(!strcmp(buf, "trw-s")) method = 0;
				if(!strcmp(buf, "bp")) method = 1;
			}
			mxFree(buf);
		}
		if((curField = mxGetField(oInPtr, 0, "maxIter")) != NULL){
			MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<maxIter>>");
			MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<maxIter>>");
			options.m_iterMax = (int)(*(double*)mxGetData(curField));
			MATLAB_ASSERT(options.m_iterMax >= 1, "Wrong value for options.maxIter: expected value is >= 1");
		}
		if((curField = mxGetField(oInPtr, 0, "verbosity")) != NULL){
			MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<verbosity>>");
			MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<verbosity>>");
			verbosityLevel = (int)(*(double*)mxGetData(curField));
			MATLAB_ASSERT(verbosityLevel == 0 || verbosityLevel == 1 || verbosityLevel == 2, "Wrong value for options.verbosity: expected value is 0, 1, or 2");
		}
		if((curField = mxGetField(oInPtr, 0, "funcEps")) != NULL){
			MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<funcEps>>");
			MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<funcEps>>");
			options.m_eps = *(double*)mxGetData(curField);
		}
		if((curField = mxGetField(oInPtr, 0, "printMinIter")) != NULL){
			MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<printMinIter>>");
			MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<printMinIter>>");
			options.m_printMinIter = (int)(*(double*)mxGetData(curField));
			MATLAB_ASSERT(options.m_printMinIter >= 0, "Wrong value for options.printMinIter: expected value is >= 0");
		}
		if((curField = mxGetField(oInPtr, 0, "printIter")) != NULL){
			MATLAB_ASSERT(mxGetClassID(curField) == mxDOUBLE_CLASS, "Wrong structure type for options: expected DOUBLE for field <<printIter>>");
			MATLAB_ASSERT(mxGetNumberOfElements(curField) == 1, "Wrong structure type for options: expected 1 number for field <<printIter>>");
			options.m_printIter = (int)(*(double*)mxGetData(curField));
			MATLAB_ASSERT(options.m_printIter >= 1, "Wrong value for options.printIter: expected value is >= 1");
		}
	}	

	// get unary potentials
	MATLAB_ASSERT(mxGetNumberOfDimensions(uInPtr) == 2, "Unary term array is not 2-dimensional");
	MATLAB_ASSERT(mxGetPi(uInPtr) == NULL, "Unary potentials should not be complex");
	
	mwSize numNodes = mxGetN(uInPtr);
	mwSize numLabels = mxGetM(uInPtr);

	MATLAB_ASSERT(numNodes >= 1, "The number of nodes is not positive");
	MATLAB_ASSERT(numLabels >= 1, "The number of labels is not positive");
	MATLAB_ASSERT(mxGetClassID(uInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for input unary term argument");
	double* termW = (double*)mxGetData(uInPtr);

	//get label matrix
	double* labelMatrix = NULL;
	if(mInPtr != NULL){
		if(!mxIsEmpty(mInPtr)){
			MATLAB_ASSERT(mxGetNumberOfDimensions(mInPtr) == 2, "Label matrix is not 2-dimensional");
			MATLAB_ASSERT(mxGetPi(mInPtr) == NULL, "Label matrix should not be complex");
			MATLAB_ASSERT(mxGetClassID(mInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for label matrix");
			MATLAB_ASSERT(mxGetN(mInPtr) == numLabels && mxGetM(mInPtr) == numLabels, "Label matrix should be of size NumLabels x NumLabels");
			
			labelMatrix = (double*)mxGetData(mInPtr);
		}
	}

	//get pairwise potentials
	MATLAB_ASSERT(mxIsSparse(pInPtr), "Expected sparse array for neighbours");
	MATLAB_ASSERT(mxGetN(pInPtr) == numNodes && mxGetM(pInPtr) == numNodes,
	              "Neighbours array must be NumNodes x NumNodes in size");
	MATLAB_ASSERT(mxGetClassID(pInPtr) == mxDOUBLE_CLASS, "Expected mxDOUBLE_CLASS for neighbours array");
	MATLAB_ASSERT(mxGetPi(pInPtr) == NULL, "Pairwise potentials should not be complex");

	mwIndex colNum = (mwIndex)mxGetN(pInPtr);
	const mwIndex* ir = mxGetIr(pInPtr);
	const mwIndex* jc = mxGetJc(pInPtr);
	double*        pr = mxGetPr(pInPtr);

	////check pairwise terms
	//mwSize numEdges = 0;
	//for (mwIndex c = 0; c < colNum; ++c) {
	//	mwIndex rowStart = jc[c]; 
	//	mwIndex rowEnd   = jc[c+1]; 
	//	for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
	//		mwIndex r = ir[ri];

	//		double dw = pr[ri];
	//		if( r < c) numEdges++;
	//	}
	//}

	
	//create MRF object
	MRFEnergy<TypeGeneral>* mrf;
	MRFEnergy<TypeGeneral>::NodeId* nodes;
	std::vector<TypeGeneral::REAL> energy, lowerBound;
	std::vector<clock_t> time_arr;
	
	TypeGeneral::REAL *D = new TypeGeneral::REAL[numLabels];
	TypeGeneral::REAL *P = new TypeGeneral::REAL[numLabels * numLabels];
	for(int i = 0; i < numLabels * numLabels; ++i)
			P[i] = 0;

	mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize());
	nodes = new MRFEnergy<TypeGeneral>::NodeId[numNodes];
	
	// construct energy
	// add unary terms
	for(int i = 0; i < numNodes; ++i){
		nodes[i] = mrf->AddNode(TypeGeneral::LocalSize(numLabels), TypeGeneral::NodeData(termW + i * numLabels));
	}

	//add pairwise terms
	for (mwIndex c = 0; c < colNum; ++c) {
		mwIndex rowStart = jc[c]; 
		mwIndex rowEnd   = jc[c + 1]; 
		for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
			mwIndex r = ir[ri];
			double dw = pr[ri];

//			if (r < c) {
//  NOT working in current TRW-S code!!!!!			
//				if (dw >= 0) 
//					mrf->AddEdge(nodes[r], nodes[c], TypeGeneral::EdgeData(TypeGeneral::POTTS, dw));
//				else
					if (labelMatrix == NULL){ //If label matrix is not specified add Potts term
						for(int i = 0; i < numLabels; ++i)
							for(int j = 0; j < numLabels; ++j)
								if(i != j) P[j + numLabels * i] = dw;
					}
					else
					{ //add matrix that is specified by user
						for(int i = 0; i < numLabels; ++i)
							for(int j = 0; j < numLabels; ++j)
								P[j + numLabels * i] = dw * labelMatrix[j + numLabels * i];
					}

					mrf->AddEdge(nodes[r], nodes[c], TypeGeneral::EdgeData(TypeGeneral::GENERAL, P));
//				}
		}
	}	

	/////////////////////// TRW-S algorithm //////////////////////
	if (verbosityLevel < 2)
		options.m_printMinIter = options.m_iterMax + 2;

	clock_t tStart = clock();
		
	if(method == 0) //TRW-S
	{
		// Function below is optional - it may help if, for example, nodes are added in a random order
		mrf->SetAutomaticOrdering();

		mrf->Minimize_TRW_S(options, lowerBound, energy, time_arr);

		if(verbosityLevel >= 1)
			printf("TRW-S finished. Time: %f\n", (clock() - tStart) * 1.0 / CLOCKS_PER_SEC);
	}
	else
	{
		// Function below is optional - it may help if, for example, nodes are added in a random order
		//mrf->SetAutomaticOrdering();

		mrf->Minimize_BP(options, energy, time_arr);
		// lowerBound = std::numeric_limits<double>::signaling_NaN();

		if(verbosityLevel >= 1)
			printf("BP finished. Time: %f\n", (clock() - tStart) * 1.0 / CLOCKS_PER_SEC);
	}

	//output the energy values array
	if(eOutPtr != NULL)	{
		*eOutPtr = mxCreateNumericMatrix(energy.size(), 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData(*eOutPtr);
		for(int i = 0; i < energy.size(); ++i)
			segment[i] = (double)(energy[i]);
	}

	//output the best solution
	if(sOutPtr != NULL)	{
		*sOutPtr = mxCreateNumericMatrix(numNodes, 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData(*sOutPtr);
		for(int i = 0; i < numNodes; ++i)
			segment[i] = (double)(mrf -> GetSolution(nodes[i])) + 1;
	}

	//output the lower bound array
	if(lbOutPtr != NULL)	{
		*lbOutPtr = mxCreateNumericMatrix(lowerBound.size(), 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData(*lbOutPtr);
		for(int i = 0; i < lowerBound.size(); ++i)
			segment[i] = (double)(lowerBound[i]);
	}

	//output the lower bound array
	if(tOutPtr != NULL)	{
		*tOutPtr = mxCreateNumericMatrix(time_arr.size(), 1, mxDOUBLE_CLASS, mxREAL);
		double* segment = (double*)mxGetData(*tOutPtr);
		for(int i = 0; i < time_arr.size(); ++i)
			segment[i] = (double)(time_arr[i]);
	}

	// done
	delete [] nodes;
	delete mrf;
	delete [] D;
	delete [] P;
}

