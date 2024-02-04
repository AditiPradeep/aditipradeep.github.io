---
layout: page
title: Select coding examples from SuperCDMS
---

I am an active software developer for SuperCDMS's data processing packages and data acquisition packages. SuperCDMS uses C++ for its core packages which handles large amounts of data for computational efficiency. SuperCDMS analyses are done in both python (numpy, pandas, scipy, matplotlib) and  C++ ([ROOT's RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) which provides a nice python interface). All SuperCDMS softwares are privately owned and contain sensitive information which cannot be fully publically released. Therefore, I'm trying to provide a couple of snippets of my ``software portfolio" to demonstrate my extensive experience coding for SuperCDMS for the past 5 years.

## Data processing packages

SuperCDMS data is acquired in [Midas format](https://daq00.triumf.ca/MidasWiki/index.php/Main_Page), which is a commonly used data format in the particle physics community. The data acquired is transferred to the [scientific High Performance Comupting facility](https://sdf.slac.stanford.edu/public/doc/#/) at the Stanford Linear Accelerator Center (SLAC) for storage and processing. SuperCDMS data processing is handled by a group of packages and builds heavily on Object Oriented Programming at its core. Each package/sub-package has a well defined role which I demonstrate below with the cartoon I made for the documentation effort I led for these packages using gitlab CI/CD. Note that I also currently maintain the Gitlab documentation.

![SuperCDMS processing packages]({{'/assets/img/Processing_code_structure.png' | relative_url}})

BatCommon is the package where data analyses classes, data reading classes and utility functions live. CDMSBats is the main package which inherits classes from BatCommon and applies it to its main processing sub-packages which I have extensively contributed to. BatFaker is a simulations package and BatViewer is a display package to view data. The constantly utilized packages are BatNoise, BatRoot and BatCalib. I led the effort to perform critical updates to these three sub-packages and BatCommon to prepare the processing software to process data when the new experiment begins data taking in 2025. We have also implemented unit testing for CDMSBats with gitlab CI/CD. Below is the rough work flow of the processing packages (note that BatCommon supplies all the necessary analyses classes and tools to CDMSBats behind the scenes).

![CDMSBats flow chart]({{'/assets/img/CDMSBats_flowchart.png' | relative_url}})

Here I describe one avenue where I contributed extensively (about 2000+ lines of code additions). The new data for SuperCDMS will be taken with two different sampling frequencies in the same waveform (check figure below). In the figure, the waveform on the left has its tails on the left and right of the peak downsampled by a factor of 16 (by averaging every 16 data points). This "hybrid" waveform you get after downsampling is shown on the right. All our regression is performed in frequency domain with Fourier transforms and this hybrid sampling business makes that difficult. This calls for a complete do-over of our regression algorithms and several other important processing algorithms. I developed the idea for the regression algorithm (called 2-speed Optimal Filter) mathematically, performed preliminary testing, implemented it into code and validated it. In addition to this algorithm, most other algorithms had to be modified to accommodate this new hybrid sampling procedure.

![From normal pulse to hybrid pulse]({{'/assets/img/Hybrid_pulse_demo.png' | relative_url}})

In the code I show below, "OP" or "On-Pulse" or "fast" refers to the non-downsampled part just around the peak. The "DS" or "Downsampled" pulse refers to what you would get if you uniformly downsampled the whole trace without leaving out the OP region. The algorithm I developed is callled the 2-Speed Optimal Filter. Here is a copy of the BatCommon class I wrote for this: [Example 2-speed Optimal Filter Class]({{ '/OptimalFilterPhononTwoSpeed.cxx' | relative_url}}). This class is used in the following CDMSBats snippet under another class called EventBuilder in BatRoot:

```C++
#include "OptimalFilterPhononTwoSpeed.h"

//Optimal filter for hybrid data
//Two-speed Optimal Filter for phonon
void EventBuilder::DoOptimalFilterPhononTwoSpeed(int detNum, const string& sensorType)
{
	.
	.
	.
	.
	.
	.
	.
//----------------------------------------
//Skipping a bunch of lines at the start of this function which loads necessary parameters such as pulses,
// BatNoise outputs from filter file (refer flow chart) and user settings.
// Jumping straight to the BatCommon class call:
//-----------------------------------------

	  // ------- set OptimalFilterPhonon  parameters ---------
	  //Initialize object
	  OptimalFilterPhononTwoSpeed tempTwoSpeedOptimalFilter; 
	  tempTwoSpeedOptimalFilter.SetWindows(win1,win2);
	  tempTwoSpeedOptimalFilter.SetOPWindowExtension(extend_OPwindowLeft_by,extend_OPwindowRight_by);
	  tempTwoSpeedOptimalFilter.SetTemplateOffset(template_offset);
	  tempTwoSpeedOptimalFilter.SetTriggerOffset(preTrigger);
	  tempTwoSpeedOptimalFilter.SetDownSamplingFactor(dsFactor);
	  tempTwoSpeedOptimalFilter.LoadOnPulseTemplates(templateFFTOP,optimalFilterConjOP, numTemplates, num_on_pulse);
	  tempTwoSpeedOptimalFilter.LoadOnPulseNormalizations(normFFTOP, sigToNoiseSqOP, noiseFFTsqOP);
	  tempTwoSpeedOptimalFilter.LoadInvCovMat(invCovMat);
	  //tempTwoSpeedOptimalFilter.LoadCutoffFreq(cutoffFreq);
	   
	  //set timing parameters
	  tempTwoSpeedOptimalFilter.SetOnPulseSampleTime(1.0/sampleRate);
	  tempTwoSpeedOptimalFilter.SetDownsampledSampleTime(1.0/sampleRateLow);

	  // -------  do Optimal Filter ------------ 
	  tempTwoSpeedOptimalFilter.DoOFWithDelayScan(opPulse);
	  tempTwoSpeedOptimalFilter.LoadDownsampledTemplates(templateTime, noiseFFTsqDS, num_pre_pulse*dsFactor);//must be done only after OnPulse optimal filtering
	  tempTwoSpeedOptimalFilter.DoCalcForDownsampled(dsPulse);
	  tempTwoSpeedOptimalFilter.CombineAmplitudes();
	  tempTwoSpeedOptimalFilter.LoadDownsampledQuantitiesForAmp0(templateFFTDS, optimalFilterConjDS, normFFTDS, sigToNoiseSqDS);
	  tempTwoSpeedOptimalFilter.DoCalcForAmp0(dsPulse);

	  //cleanup and free memory
	  delete [] templateFFTOP;
	  delete [] optimalFilterConjOP;
 
	  // -------  store  information -----------
      
	  //store instance of this class so RQs can be read out a little later
	  aPulseData->StorePulseAnalysis(tempTwoSpeedOptimalFilter);
      
      // I've not shown the beginning of these loops in this snippet
	}  // ======= end loop ZIP pulses  ======= 

    } //end if detNum found in map
  return; 
}

```
The object calls in the above snippet follow this flow chart: 

![Code snippet flow chart]({{'/assets/img/2SOF_flowchart.png' | relative_url}})

## Data Acquisition software

The SuperCDMS Data Acquisition Software (SuperCDMS DAQ) also consists of many sub-packages with various functions. Unlike the processing software, there is no strict heirarchy in these packages.

![SuperCDMS DAQ flow chart]({{'/assets/img/SuperCDMS_DAQ_flowchart.png' | relative_url}})

I am a developer for MidasDAQ and the official release manager for the umberella package SuperCDMS DAQ. I organize and coordinate timely release of the DAQ package, collecting tages for the sub-packages and performing a gitflow based version tagged release. We use the jira issue tracking system for bug reports and tracking. 

Below is an example of a stability monitoring function I wrote for the MidasDAQ. The function gets called every few seconds to query the High voltage supply which applies a voltage on the detector, checks if the current across the detector is below a preset threshold, issues an error message if the current is too large and follows up by turning off the power supply in case of emergencies. Note that a current across the detector can be quite disastrous for an expensive equipment.

![Example function]({{'/assets/img/DAQ_snippet1.png' | relative_url}})

Here is another function which retrieves the CPLD hardware version from our Detector Control and Readout Card (DCRC) and writes it to our web interface called the ODB interface. Retrieving the versions requires reading and writing to specific registers on the DCRC which you can see in the code:

![Example function]({{'/assets/img/DAQ_snippet2.png' | relative_url}})