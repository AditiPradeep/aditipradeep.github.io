////////////////////////////////////////////////////////////////////////////////                                                                                                                           
//Class Name: OptimalFilterPhononTwoSpeed  
//Author:  A. Pradeep
                                                                                                                                                                           
//Description: This class perfoms an optimal filtering using noise fft and signal fft
// on on-pulse and downsampled data from a hybrid pulse. 
//It is adapted from the original phonon optimal filtering algorithm.                    
// Original author of the darkpipe code:  R. Schnee                                                                                                                                             

//File Import By: A. Pradeep                                                                                                                                                                 

//Creation Date: Feb. 11, 2021                                                                                                                                                                             

//Modifications: 
////////////////////////////////////////////////////////////////////////////////// 

#include <iostream>
#include <limits>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TPrincipal.h"

#include "OptimalFilterPhononTwoSpeed.h"
#include "PulseTools.h"

////////////////////////////////////////////////////////                                                                                                                                                    

//do not modify the signature of this constructor                                                                                                                                                           
//instead use InitializeParameters() to pass in values to your class  
OptimalFilterPhononTwoSpeed::OptimalFilterPhononTwoSpeed(const string& className) :
  fOPAmp(-999999.),
  fOPAmp0(-999999.),
  fDSAmp0(-999999),
  fDSAmp(-999999.),
  fOPChisq(-999999.),
  fDSChisq(-999999.),
  fOPDelay(-999999.),
  fDSDelay(-999999.),
  fCombinedAmp(-999999),
  fCombinedAmp0(-999999),
  fOnPulsedT(1.6e-6),
  fDownsampleddT(25.6e-6),
  fwindow1(-999999),
  fwindow2(-999999),
  //fCutoffFreq(-999999.),
  fOnPulseTemplatesLoaded(false),
  fOnPulseNormalizationsLoaded(false),
  fDownsampledTemplatesLoaded(false),
  fDownsampledNormalizationsLoaded(false),
  fCalcCov(false),
  fOnPulseNbins(-999999),
  fDownsampledNbins(-999999),
  fOnPulseWindowExtensionLeft(-999999),
  fOnPulseWindowExtensionRight(-999999), 
  fNumOfPrePulseBins(-999999),    
  fTemplateOffset(-999999),
  fTriggerOffset(-999999),
  fDownSamplingFactor(-999999)
{

  //   cout <<"Hello from OptimalFilterPhonon()" << endl;                                                                                                                                                  

  //these members along with fRQlist are inherited from TCDMSAnalysis                                                                                                                                      
  fClassName = className;
  fStoreRQs = true;
  /*
   if (fClassName.find("OptimalFilterPhononGlitch1")!=string::npos)
        fAnalysisInitials = "glitch1OF";
   else if (fClassName.find("OptimalFilterPhononPileup")!=string::npos)
        fAnalysisInitials = "pileupOF";
   else if (fClassName.find("OptimalFilterPhononLFnoise1")!=string::npos)
        fAnalysisInitials = "lfnoise1OF";
   else if (fClassName.find("OptimalFilterPhononDMC")!=string::npos)
        fAnalysisInitials = "dmcOF";
	else*/
  fAnalysisInitials = "TwoSOF";


  //Construct the RQ list                                                                                                                                                                                  
  ConstructRQList();

  //initialization of member data can go here                                                              

}

OptimalFilterPhononTwoSpeed::~OptimalFilterPhononTwoSpeed()
{
  //   cout <<"Goodbye from OptimalFilterPhonon()" << endl;                                             
}
//This method constructs the RQ list that is handed off to BatOutputManager  
//It also sets the default value for the RQ to initVal.                                                                                                                                                   

void OptimalFilterPhononTwoSpeed::ConstructRQList()
{
  double initVal = -999999.;

  //construct the RQ list here                                                                                                                                                                            
  //OP=On-Pulse, DS=DownSampled
  
   fRQList.insert(pair<string,double>( fAnalysisInitials + "OPAmp",  initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "DSAmp", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "OPChisq", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "DSChisq", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "OPDelay", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "DSDelay", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "CombinedAmp", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "OPAmp0", initVal));   
   fRQList.insert(pair<string,double>( fAnalysisInitials + "DSAmp0", initVal));
   fRQList.insert(pair<string,double>( fAnalysisInitials + "CombinedAmp0", initVal));
  //Any RQ that is included in the above list will be written out by BatRoot.  Add to this as you please.

  return;
}
//This is the main call
 
void OptimalFilterPhononTwoSpeed::DoOFWithDelayScan(const vector<double>& aPulse)
{
  //do your calculation here!                                                                                                                                                                             
  //for debugging only                                                                                                                                                                                    
  //ConstructFakePulse(0.5e9, 525);                                                                                                                                                                       
  //This is general purpose optimal filter that scans for the minimum chi square (instead of maximum amplitude).
  //The code writes the amplitude and delay corresponding to minimum chi square as RQs.
  //The scanning is done manually instead of by an inverse Fourier transform.

  int nBins = aPulse.size();
  double sqrtdT = sqrt(fOnPulsedT); //for efficiency                                                                                                                                                      

  //============== Some Preliminary Checks =================                                                                                                                                              
  if(!fOnPulseTemplatesLoaded)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::DoOFWithDelayScan ERROR!  attempting to run OptimalFilter without loading templates!" << endl;
      exit(1);
    }

  //checking that the lengths of templates agree with the pulse lengths                                                                                                                                   
  if(nBins != fOnPulseNbins)
    {
      cerr<<"OptimalFilterPhononTwoSpeed::ERROR! Number of bins in Onpulse does not match number of bins in templates!" << endl;
      exit(1);
    }

  //================== Calculations =====================                                                                                                                                                 
  TComplex comp_zero(0.,0.);

  //vectors to hold intermediate calculations                                                                                                                                                          
  vector<TComplex> pulseFFT;

  //PulseTools::RealToComplexFFT(fFakePulse, re_PulseFFT, im_PulseFFT); //for testing only                               
  //===== 1. Calculate amplitude for each template =====
  PulseTools::RealToComplexFFT(aPulse, pulseFFT);
  
  vector<double> scannedAmp;
  vector<double> scannedChisq ;
  int numTemplates = fOnPulseTemplateFFT.size();
  
  for (int templItr = 0; templItr < numTemplates; templItr++)
    {
      vector<TComplex> pProd;
      double sum_amp = 0.0;
      for(int binItr=0; binItr < nBins; binItr++)
	{
	  if (templItr == 0)
	    {
	      pulseFFT[binItr] *= sqrtdT;
	    }
	  pProd.push_back(pulseFFT[binItr]*fOnPulseOptimalFilter[templItr][binItr]);
	     
	  if(binItr != 0) 
	    {
	      sum_amp += pProd[binItr].Re();
	    }
	  //if(templItr == numTemplates-1){cout<<sum_amp<<endl;}
	}
  
      scannedAmp.push_back(sum_amp/fOnPulseSigToNoiseSq[templItr]);
      pProd.clear();
    }
  //======== 2. Calculate chisq for each template =========-
  for (int templItr= 0 ; templItr < numTemplates; templItr++)
    {
      double chisq=0;
      for (int binItr = 1; binItr < nBins; binItr++){
	TComplex fit_fft( scannedAmp[templItr]*fOnPulseTemplateFFT[templItr][binItr] );
  
	double chisqBin =  pow(TComplex::Abs(pulseFFT[binItr] - fit_fft), 2) / fOnPulseNoiseFFTSq[binItr];
	chisq += chisqBin;
      }

      scannedChisq.push_back(chisq);
    }
  /*
   for (int i =0; i<numTemplates; i++)
   { cout<<"Scanned amp"<<scannedAmp[i]<<endl;
   cout<<"Scanned chi"<<scannedChisq[i]<<endl;
     }
  */
  //===== 3. get min chisq  and delay =====                                                                                                                                                             

  double minChisq = 1*numeric_limits<double>::infinity();     
  double finalAmp = 0.0;
  int idelay = 0;

  //pick out minimum chisq and corresponding delay 
  for(int templItr = 0; templItr < numTemplates; templItr++) 
    {
      if(scannedChisq[templItr] < minChisq) 
	{
	  minChisq = scannedChisq[templItr];
	  finalAmp = scannedAmp[templItr]; 
	  idelay = templItr - fOnPulseWindowExtensionLeft; 
	}
    }
  //================== Store Results ====================                                                                                                                                                 
  fOPAmp = finalAmp;
  fOPAmp0 = scannedAmp[fOnPulseWindowExtensionLeft];
  fOPChisq = minChisq;
  fOPDelay = idelay*fOnPulsedT;
  // cout<<"OPAmp="<<fOPAmp<<endl;
  //cout<<"OPchi="<<fOnPulseChisq<<endl;
  //cout<<"OPdelay="<<fOPDelay<<endl;
  //Next, store the results of this calculation as the RQ's.                                                                                                                                              
  //These values will be included in the output of BatRoot.                                                                                                                                               
  if(fStoreRQs && !fCalcCov) {
    fRQList[fAnalysisInitials + "OPAmp"] = fOPAmp;
    fRQList[fAnalysisInitials + "OPAmp0"] = fOPAmp0;
    fRQList[fAnalysisInitials + "OPChisq"] = fOPChisq;
    fRQList[fAnalysisInitials + "OPDelay"] = fOPDelay;
  }

  //============= Delete Templates to minimize copying  ================                                                                                                                                  
  if(!fCalcCov)
    {
      //fOnPulseTemplate.clear();
      fOnPulseNoiseFFTSq.clear();
      fOnPulseTemplateFFT.clear();
      fOnPulseOptimalFilter.clear();
       
      // ========== cleanup! =========== 
      fOnPulseTemplatesLoaded = false;
      fOnPulseNormalizationsLoaded = false;
    }
  
  return;
}

void OptimalFilterPhononTwoSpeed::DoCalcForDownsampled(const vector<double>& aPulse)
{
  //do your calculation here!                                                                                                                                                                              

  //for debugging only                                                                                                                                                                                     
  //ConstructFakePulse(0.5e9, 525);                                           

  //Normal optimal filter code from OptimalFilterPhonon.cxx                                                                                                                             

  int nBins = aPulse.size();
  double sqrtdT = sqrt(fDownsampleddT); //for efficiency

  //============== Some Preliminary Checks =================                                                                                                                                               

  if(!fDownsampledTemplatesLoaded)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::DoCaclForDownsampled ERROR!  attempting to run OptimalFilter without loading templates!" << endl;
      exit(1);
    }

  //checking that the lengths of templates agree with the pulse lengths                                                                                                                                   
  if(nBins != fDownsampledNbins)
    {
      cerr<<"OptimalFilterPhononTwoSpeed::ERROR! Number of bins in Downsampledpulse does not match number of bins in templates!" << endl;
      exit(1);
    }

  //checking that the fit window is valid                                                                                                                                                                  
  if(fwindow1 < 0 || fwindow2 < 0)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::ERROR! Fit window appears to be uninitialized!" << endl;
      exit(1);
    }
  //================== Calculations =====================                                                                                                                                                  
  TComplex comp_zero(0.,0.);
  //vectors to hold intermediate calculations                                                                                                                                                              
  vector<TComplex> pProd;
  vector<TComplex> pulseFFT;
  vector<double> p_prod_ifftRe;

  //===== 1. construct fitpulse fft =====                                                                                                                                                                  
  //PulseTools::RealToComplexFFT(fFakePulse, re_PulseFFT, im_PulseFFT); //for testing only                                                                                                                 
  PulseTools::RealToComplexFFT(aPulse, pulseFFT);

  double amp0 = 0.;
  //===== 2. construct p_prod = sqrt(fdT) * pulse_fft * s*/J =====
  for(int binItr=0; binItr < nBins; binItr++)
    {
      pulseFFT[binItr] *= sqrtdT;
      pProd.push_back(pulseFFT[binItr]*fDownsampledOptimalFilter[binItr]);

      if(binItr != 0) 
	{amp0 += pProd[binItr].Re();}
    }

  //normalizing amp0's                                                                                                                                                                                     
  amp0 /= fDownsampledSigToNoiseSq;
  //cout<<"Amp0="<<amp0<<endl;
  //ignoring DC component;                                                                                                                                                                                 
  pProd[0] = comp_zero;
  //===== 3. get max amplitude and delay =====                                                                                                                                                             
  
  vector<double> p_ahat; //vector of all amplitudes                                                                                                                                                        
  double maxAmp = -1*numeric_limits<double>::infinity(); //maxAmp can be a negative number                                                                                                                 
  double finalAmp = 0.0;
  int idelay = 0;

  PulseTools::ComplexToRealIFFT(pProd, p_prod_ifftRe);  //C2R - no imaginary component!                                                                                                                    

  //pick out maximum amplitude and corresponding delay                                                                                                                                                     
  for(int binItr=0; binItr < nBins; binItr++) {

    double amp =  p_prod_ifftRe[binItr];

    if(amp > maxAmp) {
      maxAmp = amp;
      finalAmp = p_prod_ifftRe[binItr]/fDownsampledNormFFT;
      idelay = binItr; }
    //search up to fwindow1, and after fwindow2                                                                                                                                                            
    if(binItr == fwindow1 - 1)
      binItr = fwindow2 - 2; //-2 b/c for loop increments by 1 before next round, and c++ array convention adds -1                                                                                       
  }
  
  //compute phase factor with this delay                                                                                                                                                                  
  int delay = (idelay < nBins/2 ? idelay : (idelay - nBins));
  //cout<<delay<<endl;
  //cout<<"finalAmp="<<finalAmp<<endl;
  //cout<<"Norm="<<fDownsampledNormFFT<<endl;
  //===== 4. compute chisq =====                                                                                                                                                                           

  double chisq = 0;

  //ignoring DC component                                                                                                                                                                                  
  for(int binItr=1; binItr < nBins; binItr++)
    {
      double theta = 2.0*TMath::Pi()*((double)binItr/(double)nBins)*(double)delay;

      TComplex phase_factor(cos(theta), sin(theta));
      TComplex fit_fft( finalAmp*fDownsampledTemplateFFT[binItr]/phase_factor );

      double chisqBin =  pow(TComplex::Abs(pulseFFT[binItr] - fit_fft), 2)/fDownsampledNoiseFFTSq[binItr];
      chisq += chisqBin;
    }
  //================== Store Results ====================                                                                                                                                                  

  //this is a little redundant, do I want to keep this?                                                                                                                                                    
  fDSAmp = finalAmp;
  fDSChisq = chisq;
  fDSDelay = delay*fDownsampleddT;
  //cout<<"DSAmp="<<fDSAmp<<endl;
  //cout<<"DSdelay=" <<fDSDelay<<endl;
  //cout<<"DSchi="<<fDSChisq<<endl;
  //Next, store the results of this calculation as the RQ's.                                                                                                                                               
  //These values will be included in the output of BatRoot.                                                                                                                                                
  if(fStoreRQs && !fCalcCov) {
    fRQList[fAnalysisInitials + "DSAmp"] = fDSAmp;
    fRQList[fAnalysisInitials + "DSChisq"] = fDSChisq;
    fRQList[fAnalysisInitials + "DSDelay"] = fDSDelay;
  }

  //============= Delete Templates to minimize copying  ================                                                                                                                                   

  //fPulseTemplate.clear();
  fDownsampledNoiseFFTSq.clear();
  fDownsampledTemplateFFT.clear();
  fDownsampledOptimalFilter.clear();
  
  // ========== cleanup! ===========                                                                                                                                                                       
  fDownsampledTemplatesLoaded = false;
  fDownsampledNormalizationsLoaded = false;
  //cout<<"Exiting function"<<endl;
  return;
}
void OptimalFilterPhononTwoSpeed::DoCalcForAmp0(const vector<double>& aPulse)
{
  //do your calculation here!                                                                                                                                                                              

  //for debugging only                                                                                                                                                                                     
  //ConstructFakePulse(0.5e9, 525);                                           

  //Normal optimal filter code from OptimalFilterPhonon.cxx                                                                                                                             

  int nBins = aPulse.size();
  double sqrtdT = sqrt(fDownsampleddT); //for efficiency

  //============== Some Preliminary Checks =================                                                                                                                                               

  if(!fDownsampledQuantities0Loaded)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::DoCaclForAmp0 ERROR!  attempting to run OptimalFilter without loading templates!" << endl;
      exit(1);
    }

  //checking that the lengths of templates agree with the pulse lengths                                                                                                                                   
  if(nBins != fDownsampledNbins)
    {
      cerr<<"OptimalFilterPhononTwoSpeed::ERROR! Number of bins in Downsampledpulse does not match number of bins in templates!" << endl;
      exit(1);
    }

  //================== Calculations =====================                                                                                                                                                  
  TComplex comp_zero(0.,0.);
  //vectors to hold intermediate calculations                                                                                                                                                              
  vector<TComplex> pProd;
  vector<TComplex> pulseFFT;
  vector<double> p_prod_ifftRe;

  //===== 1. construct fitpulse fft =====                                                                                                                                                                  
  //PulseTools::RealToComplexFFT(fFakePulse, re_PulseFFT, im_PulseFFT); //for testing only                                                                                                                 
  PulseTools::RealToComplexFFT(aPulse, pulseFFT);
  double amp0 = 0.;
  //===== 2. construct p_prod = sqrt(fdT) * pulse_fft * s*/J =====
  for(int binItr=0; binItr < nBins; binItr++)
    {
      pulseFFT[binItr] *= sqrtdT;
      pProd.push_back(pulseFFT[binItr]*fDownsampledOptimalFilter0[binItr]);

      if(binItr != 0) 
	{amp0 += pProd[binItr].Re();}
    }

  //normalizing amp0's                                                                                                                                                                                     
  amp0 /= fDownsampledSigToNoiseSq0;
  fDSAmp0 = amp0;
  //================= Calculate combinedamp0 =============

  double numerator = (fCovMatInv[0][0] + fCovMatInv[0][1])*fOPAmp0 + (fCovMatInv[1][0] + fCovMatInv[1][1])* fDSAmp0;
  double denominator = (fCovMatInv[0][0] + 2*fCovMatInv[0][1] + fCovMatInv[1][1]);

  fCombinedAmp0 = numerator/denominator;

  //================== Store Results ====================                                                                                                                                                  

  //cout<<"DSAmp="<<fDownsampledAmp<<endl;
  //cout<<"DSdelay=" <<delay*fDownsampleddT <<endl;
  //cout<<"DSchi="<<fDownsampledChisq<<endl;
  //Next, store the results of this calculation as the RQ's.                                                                                                                                               
  //These values will be included in the output of BatRoot.                                                                                                                                                
  if(fStoreRQs && !fCalcCov) {
    fRQList[fAnalysisInitials + "DSAmp0"] = fDSAmp0;
    fRQList[fAnalysisInitials + "CombinedAmp0"] = fCombinedAmp0;
  }

  //============= Delete Templates to minimize copying  ================                                                                                                                                   

  //fPulseTemplate.clear();
  fDownsampledTemplateFFT0.clear();
  fDownsampledOptimalFilter0.clear();
 
  // ========== cleanup! ===========                                                                                                                                                                       
  fDownsampledQuantities0Loaded = false;
  //cout<<"Exiting function"<<endl;
  return;
}

void OptimalFilterPhononTwoSpeed::LoadOnPulseTemplates(const TComplex* pulseTemplateFFT, const TComplex* optimalFilter, int numTemplates, int opSize)
{
  //these templates should all have the same lengths                                                                                                                                                   
 
  fOnPulseNbins = opSize;
  /*
   if((int)optimalFilter[0].size() != fOnPulseNbins)
     {
     cerr <<"OptimalFilterPhononTwoSpeed::ERROR!  Template lengths do not match, check the input to LoadOnPulseTemplates." << endl;
       exit(1);
       }*/

  if(fOnPulseWindowExtensionLeft+fOnPulseWindowExtensionRight+1 != numTemplates)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::ERROR!  Mismatch in number of templates in, check the input to LoadOnPulseTemplates." << endl;
      exit(1);
    }

  for (int itr = 0; itr < numTemplates; itr ++)
    {
      vector<TComplex> tempTemplateVec;
      vector<TComplex> tempOFVec;
      for (int binItr = itr*opSize; binItr < (itr+1)*opSize; binItr++)
	{
	  tempTemplateVec.push_back(pulseTemplateFFT[binItr]);
	  tempOFVec.push_back(optimalFilter[binItr]);
	}
       
      fOnPulseTemplateFFT.push_back(tempTemplateVec);
      fOnPulseOptimalFilter.push_back(tempOFVec);
      tempTemplateVec.clear();
      tempOFVec.clear();
    }

  fOnPulseTemplatesLoaded = true;

  return;
}

void OptimalFilterPhononTwoSpeed::LoadDownsampledTemplates(const vector<double>& pulseTemplateTime, const vector<double>& noiseFFTsq, int num_pre)
{
  //these templates should all have the same lengths                                                                                                                                                       
  if(fDownsampledTemplatesLoaded)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::LoadDownsampledTemplates ERROR!  Templates already loaded!" << endl;
      exit(1);
    }

  fDownsampledNbins = (int)pulseTemplateTime.size()/fDownSamplingFactor;
  fDownsampledNoiseFFTSq = noiseFFTsq;
  fNumOfPrePulseBins = (int)num_pre/fDownSamplingFactor;
  /*
   if((int)optimalFilter.size() != fOnPulseNbins)
     {
     cerr <<"OptimalFilterPhononTwoSpeed::ERROR!  Template lengths do not match, check the input to LoadOnPulseTemplates." << endl;
       exit(1);
     }
  */
  // Note to self: Positive delay means signal is ahead of the template, negative means signal is lagging behind template
  vector<double> tempDownsampledTemplateTime;
  double delayBin = (fOPDelay / fOnPulsedT);
  int PulseStartBin = num_pre + fTriggerOffset;
  
  vector<double> tempTemplateTime;
  
  // Shift the template to align with the data
  for (int binItr = 0; binItr < (int)pulseTemplateTime.size(); binItr++)
    {
      if ((binItr < PulseStartBin - delayBin) || binItr >=(int)pulseTemplateTime.size() - delayBin)  
	{tempTemplateTime.push_back(0.0);}
       
      else
	{tempTemplateTime.push_back(pulseTemplateTime[binItr + delayBin]);}
    } 
  //Downsample shifted template
  tempDownsampledTemplateTime = PulseTools::DownsamplePulse(tempTemplateTime, fDownSamplingFactor);
  fDownsampledTemplateTime = tempDownsampledTemplateTime;
  
  //tempDownsampledTemplateTime = PulseTools::Normalize(tempDownsampledTemplateTime, PulseTools::MaxADC(tempDownsampledTemplateTime));//Renormalize template
  //Compute and store FFT 
  PulseTools::RealToComplexFFT(tempDownsampledTemplateTime, fDownsampledTemplateFFT);
  
  for (uint binCtr=0; binCtr < fDownsampledTemplateFFT.size(); binCtr++)
    {
      TComplex scalefactor(sqrt(fDownsampleddT),0.0);
      fDownsampledTemplateFFT[binCtr] *= scalefactor;
    }
  //cout<<"TemplFFTScaleFactor="<<fDownsampledTemplateFFT[200]<<endl;
  vector<TComplex> tempDownsampledOF;
  TComplex dc(0.0, 0.0);
  tempDownsampledOF.push_back(dc);

  for(uint binCtr = 1; binCtr < fDownsampledNoiseFFTSq.size(); binCtr++){
    tempDownsampledOF.push_back(TComplex::Conjugate(fDownsampledTemplateFFT[binCtr])/fDownsampledNoiseFFTSq[binCtr]);}

  fDownsampledOptimalFilter = tempDownsampledOF;

  // === compute and store normalization for amplitude estimator ===                                                                                                                                      
  double tempDownsampledNormFFT = 0.0;

  //compute sum, ignoring the DC component                                                                                                                                                                
  for(uint binCtr=1; binCtr < fDownsampledNoiseFFTSq.size(); binCtr++){
    tempDownsampledNormFFT += pow(TComplex::Abs(fDownsampledTemplateFFT[binCtr]),2)/fDownsampledNoiseFFTSq[binCtr];}

  // Note, normFFT and sigToNoiseSq are really the same quantity except for 1/sqrt(N) factor 
  // we're storing both for legacy reasons now (they must be matched with terms that have  
  // the correct sqrt(N) factors in the numerator).  This can be cleaned up in the future.                                                                                                      
  fDownsampledSigToNoiseSq = tempDownsampledNormFFT;
  //cout<<"Sig2Noise="<<fDownsampledSigToNoiseSq<<endl;
  tempDownsampledNormFFT *= 1.0/sqrt(fDownsampledNoiseFFTSq.size()); // scale by 1/sqrt(N) for normFFT                                                                                                   
  fDownsampledNormFFT = tempDownsampledNormFFT;
  fDownsampledTemplatesLoaded = true;
  fDownsampledNormalizationsLoaded = true;

  tempDownsampledTemplateTime.clear();
  tempTemplateTime.clear();
  return;
}


void OptimalFilterPhononTwoSpeed::LoadOnPulseNormalizations(const vector<double>& normFFT, const vector<double>& sigToNoiseSq, const vector<double>& noiseFFTSq)
{
  //check that vector lengths match
  if((int)noiseFFTSq.size() != fOnPulseNbins)
    {
      cerr <<"OptimalFilterPhononTwoSpeed::ERROR!  Template lengths do not match, check the input to LoadOnPulseNormalizations." << endl;
      exit(1);
    }


  fOnPulseNormFFT = normFFT;
  fOnPulseSigToNoiseSq = sigToNoiseSq;
  fOnPulseNoiseFFTSq = noiseFFTSq;


  fOnPulseNormalizationsLoaded = true;

  return;
}

void OptimalFilterPhononTwoSpeed::LoadDownsampledQuantitiesForAmp0(const vector<TComplex>&templateFFT, const vector<TComplex>& optimalFilter, const double& normFFT, const double& sigToNoiseSq)
{
  if((int)optimalFilter.size() != (int)templateFFT.size())
    {
      cerr <<"OptimalFilterPhononTwoSpeed::ERROR!  Template lengths do not match, check the input to LoadDownsampledQuantitiesForAmp0()." << endl;
      exit(1);
    }

  fDownsampledTemplateFFT0 = templateFFT;
  fDownsampledOptimalFilter0 = optimalFilter;
  fDownsampledNormFFT0 = normFFT;
  fDownsampledSigToNoiseSq0 = sigToNoiseSq;
  
  fDownsampledQuantities0Loaded = true;

  return;
}

//=====================Code to calculate covariances and combined amplitude======================
void OptimalFilterPhononTwoSpeed::CalcCovarianceMatrix(const vector<vector<double>>& amps)
{
  //Pass amps as On-pulse in the first row and downsampled in the second row just to be able to keep track of things.
  // We will use TPrincipal class to calculate the covariance matrix. This needs the data points (amps in the definition above) in an array format.
  //Each row of the array should be a variable and each column represents a single measurement of all the variables
  //The pairs of amplitudes (one measurement) are passed one at a time to the TPrincipal class.
  /*
  //Pass values to the class
  double amp_array[2]={0.0};
  TPrincipal P(2); 
  int cols = (int)amps[0].size();
  for (int Itr = 0; Itr < cols; Itr++)
    {
      amp_array[0] = amps[0][Itr];
      amp_array[1] = amps[1][Itr];
      P.AddRow(amp_array);
    }

  //Calculate covariance matrix
  const TMatrixD *cov = P.GetCovarianceMatrix();

  //The matrix that comes out of this has only the lower triangular part stored, the upper triangle is set to zero.
  //To invert the matrix we need to make it symmetric as follows:
  TMatrixD cov1 = *cov;
  cov1.T(); //Get the transpose
  TMatrixDDiag cov_diag(cov1); //Pass the diagonal elements to cov_diag
  cov_diag = 0; // Set the diagonal elements of cov1 to zero. So now all elements except the upper triangular ones are zero.
  */
  //Now add
  //TMatrixD cov_symm = *cov + cov1;
  int cols = (int)amps[0].size();
  fCovMat[0][0]= PulseTools::RMS(amps[0])*PulseTools::RMS(amps[0])*cols;
  fCovMat[1][1]= PulseTools::RMS(amps[1])*PulseTools::RMS(amps[1])*cols;
  fCovMat[0][1]= PulseTools::Covariance(amps[0],amps[1]);
  fCovMat[1][0]= fCovMat[0][1];
  cout<<"Printing covariance matrix"<<endl;
  fCovMat.Print();
  return;

} 

void OptimalFilterPhononTwoSpeed::CalcInvCovarianceMatrix()
{
  //The cov mat has very small powers of ten for our data. Taking its inverse is difficult since the determinant comes out close to zero in order of magnitude.
  //Root always thinks the matrix is singular for this reason.
  //To deal with this, we multiply with 1e12 to make the numbers larger while taking the inverse and re-multiply this factor after the inverse to get the true inverse.
  //The principle is a^-1 = b*[a*b]^-1, where 'b' is a double and 'a' is our cov mat.
  //Calculate inverse covariance matrix
  TMatrixD* tempCov = &fCovMat; 
  tempCov->operator*=(1e12);
  //print error if matrix is singular
  if (tempCov->Determinant() == 0)
    {
      cout<< "OptimalFilterPhononTwoSpeed::WARNING! Covariance matrix passed to CalcInvCovarianceMatrix is singular" << endl;
      exit(1);
    }
  TMatrixD tempInv =tempCov->Invert();
  tempInv.operator*=(1e12);
  fCovMatInv += tempInv;
  cout<<"Printing inverse covariance matrix"<<endl;
  fCovMatInv.Print();
  return;

}

void OptimalFilterPhononTwoSpeed::LoadInvCovMat(const vector<double>& InvMatVector)
{
  fCovMatInv[0][0] = InvMatVector[0];
  fCovMatInv[0][1] = InvMatVector[1];
  fCovMatInv[1][0] = InvMatVector[2];
  fCovMatInv[1][1] = InvMatVector[3];
}

void OptimalFilterPhononTwoSpeed::CombineAmplitudes()
{
  //Amps must be written in the same row order as they were passed to calculate covariance matrix.
  //Write the Inverse cov mat to a vector for easy manipulation
  /*vector <vector<double>> tempInvCov;
  TMatrixD tempInv = fCovMat.Invert();
  for (int i = 0; i < 2; i++)
    {
      vector <double> temp_row;
      TVectorD row(2);
      row = TMatrixDRow(tempInv,i); 
      for (int j = 0; j < 2; j++)
	{
	  temp_row.push_back(row[j]);
	}
      tempInvCov.push_back(temp_row);
      temp_row.clear();
      }*/

  //Combine amplitudes using analytical formula
  double numerator = (fCovMatInv[0][0] + fCovMatInv[0][1])*fOPAmp + (fCovMatInv[1][0] + fCovMatInv[1][1])* fDSAmp;
  double denominator = (fCovMatInv[0][0] + 2*fCovMatInv[0][1] + fCovMatInv[1][1]);

  fCombinedAmp = numerator/denominator;
  //Next, store the result of this calculation as the RQ.                                                                                                                                               
  //These values will be included in the output of BatRoot.                                                                                                                                                
  if(fStoreRQs) {
    fRQList[fAnalysisInitials + "CombinedAmp"] = fCombinedAmp;
  }
  return;
}
//================= Temporary code for testing an development ===========================                                                                                                                   
/*
 void OptimalFilterPhononTwoSpeed::ConstructFakePulse(double norm, int delay)
 {
 if(!fTemplatesLoaded) { cerr <<"ERROR::Forgot to load templates!" << endl; }

   //Temp! using the template pulse as the data                                                                                                                                                           
   TFile f("testing_only/170319_1616_F0002_noisetrace.root");
if(!f.IsOpen()) { cerr <<"ERROR opening test noise file, check path!" << endl; exit(1); }
TTree* pulseTree = (TTree*)f.Get("pulseTree");

const int nBins = 2048; //for soudan data only!                                                                                                                                                          

//Set Branch Addresses                                                                                                                                                                                   
float noise[nBins];
pulseTree->SetBranchAddress("bsnPulse", noise);
pulseTree->GetEntry(0);
f.Close();

//Copy template into vector and histogram                                                                                                                                                                
for(int binItr = 0; binItr < nBins; binItr++)
  {
if(binItr < delay)
  fFakePulse.push_back((double)noise[binItr]);
 else
   fFakePulse.push_back((double)noise[binItr] + norm*fPulseTemplate[binItr-delay]);
}

return;
}*/
