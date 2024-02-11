---
layout: page
title: Select coding examples from SuperCDMS
---

I am an active software developer for SuperCDMS's data processing packages and data acquisition packages. SuperCDMS uses C++ for its core packages which handles large amounts of data for computational efficiency. SuperCDMS analyses are done in both python (numpy, pandas, scipy, matplotlib) and  C++ ([ROOT's RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) which provides a nice python interface). All SuperCDMS softwares are privately owned and contain sensitive information which cannot be fully publically released. Therefore, I'm trying to provide a couple of snippets of my "software portfolio" instead to demonstrate my extensive experience coding for SuperCDMS for the past 5 years.

# Python coding experience with SuperCDMS

Our data is noisy and to study reconstruction algorithms, one often has to simulate realistic noise. If one constructs the Power Spectral Density (PSD) of the noise, Monte Carlo sampling can be used to sample noise waveforms from the PSD. If there is correlation in noise between a few detector channels of interest a [cross PSD](https://en.wikipedia.org/wiki/Spectral_density#Cross_power_spectral_density) can be constructed and similarly sampled from to simulate more realistic noise waveforms. Here is a python package which wrote that performs this function: [NoiseGen repository](https://gitlab.com/AditiPradeep/noisegen).

Most of our data analysis is run in jupyter environments on scientific computing facilities using common python libraries and some SuperCDMS specific libraries. Common tasks include regression of various quantities, optimization and data enrichment. Here is an example of a Machine learning task I performed: [ML for dark matter searches]({{ '/dark_matter'| relative_url}})

# C++ coding experience with SuperCDMS

## Data processing packages

SuperCDMS data is acquired in [Midas format](https://daq00.triumf.ca/MidasWiki/index.php/Main_Page), which is a commonly used data format in the particle physics community. The data acquired is transferred to the [scientific High Performance Comupting facility](https://sdf.slac.stanford.edu/public/doc/#/) at the Stanford Linear Accelerator Center (SLAC) for storage and processing. SuperCDMS data processing is handled by a group of packages and builds heavily on Object Oriented Programming at its core. Each package/sub-package has a well defined role which I demonstrate below with the cartoon I made for the documentation effort I led for these packages using gitlab CI/CD. Note that I also currently maintain the Gitlab documentation.

![SuperCDMS processing packages]({{'/assets/img/Processing_code_structure.png' | relative_url}})

BatCommon is the package where data analyses classes, data reading classes and utility functions live. CDMSBats is the main package which inherits classes from BatCommon and applies it to its main processing sub-packages which I have extensively contributed to. BatFaker is a simulations package and BatViewer is a display package to view data. The commonly utilized sub-packages are BatNoise, BatRoot and BatCalib which convert the raw data into human-interpretable, physical quantities. I spearheaded the initiative to implement essential updates to three sub-packages along with BatCommon, ensuring the processing software is prepared for handling data when the new experiment commences data collection in 2025. We have also implemented unit testing for CDMSBats with gitlab CI/CD. Below is the rough work flow of the processing packages from left to right in that order (note that BatCommon supplies all the necessary analyses classes and tools to CDMSBats behind the scenes).

![CDMSBats flow chart]({{'/assets/img/CDMSBats_flowchart.png' | relative_url}})

**Example 1:**

**Here I describe one avenue where I contributed extensively (about 2000+ lines of code additions).** The new data acquisition process for SuperCDMS involves utilizing two different sampling frequencies within the same waveform, as illustrated in the diagram below. Specifically, the waveform on the left may have its tails downsampled by a factor of 16, resulting in the "hybrid" waveform depicted on the right. As our regression techniques primarily operate in the frequency domain with Fourier-transformed quantities, this hybrid sampling approach poses significant challenges, necessitating a comprehensive overhaul of our regression and processing algorithms.

I conceptualized the regression algorithm, dubbed the "2-speed Optimal Filter," devised its mathematical framework, conducted preliminary testing, implemented it in code, and validated its efficacy. Furthermore, numerous other algorithms required modification to accommodate this novel hybrid sampling methodology.

![From normal pulse to hybrid pulse]({{'/assets/img/Hybrid_pulse_demo.png' | relative_url}})


The overall algorithm follows this flowchart. In the code I show below, "OP", "On-Pulse" or "fast" refers to the non-downsampled part just around the peak. The "DS" or "Downsampled" pulse refers to what you would get if you uniformly downsampled the whole trace without leaving out the OP region. The algorithm I developed is callled the 2-Speed Optimal Filter.

![Code snippet flow chart]({{'/assets/img/2SOF_flowchart.png' | relative_url}})

Here is a copy of the BatCommon class I wrote for this: [Example 2-speed Optimal Filter Class]({{ '/OptimalFilterPhononTwoSpeed.cxx' | relative_url}}). This class is used in the following CDMSBats snippet under another class called EventBuilder in BatRoot:

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

## Data Acquisition software

The SuperCDMS Data Acquisition Software (SuperCDMS DAQ) also consists of many sub-packages with various functions. Unlike the processing software, there is no strict heirarchy in these packages.

![SuperCDMS DAQ flow chart]({{'/assets/img/SuperCDMS_DAQ_flowchart.png' | relative_url}})

I am a developer for MidasDAQ and the official release manager for the umberella package SuperCDMS DAQ. I organize and coordinate timely release of the DAQ package, collecting tages for the sub-packages and performing a gitflow based version tagged release. We use the jira issue tracking system for bug reports and tracking. 

**Example 2:**
Below is an example of a stability monitoring function I wrote for the MidasDAQ. The function gets called every few seconds to query the High voltage supply which applies a voltage on the detector, checks if the current across the detector is below a preset threshold, issues an error message if the current is too large and follows up by turning off the power supply in case of emergencies. Note that a current across the detector can be quite disastrous for an expensive equipment.

![Example function]({{'/assets/img/DAQ_snippet1.png' | relative_url}})

**Example 3**
Here is another function which retrieves the CPLD hardware version from our Detector Control and Readout Card (DCRC) and writes it to our web interface called the ODB interface. Retrieving the versions requires reading and writing to specific registers on the DCRC which you can see in the code:

![Example function]({{'/assets/img/DAQ_snippet2.png' | relative_url}})