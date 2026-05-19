'''Script to create histograms for FONLL cross-sections used to generate
kinematics. Cross-sections are obtained from
https://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html and this script
needs two of them, as a pT and eta which are inputs and script outputs root
file with necessary histograms. Details about individual settings are in
separate README.md in corresponding directory.'''

import argparse
import numpy
import yaml
import ROOT

parser = argparse.ArgumentParser( prog='FONLL histogram maker',
                    description='' )

parser.add_argument('-pt')
parser.add_argument('-eta')
parser.add_argument('-out')

args = parser.parse_args()

print(args)

datapt = numpy.loadtxt( args.pt )
dataeta = numpy.loadtxt( args.eta )

outFile = ROOT.TFile.Open( args.out, 'RECREATE' )

hheta = ROOT.TH1D( 'eta', '', 1200, -8, 8 )
hhpt = ROOT.TH1D( 'pT', '', 1200, 0, 300 )

for ii in range(1, 1200):
  xbin = hheta.FindBin(dataeta[ii-1][0])
  hheta.SetBinContent(xbin, dataeta[ii-1][1])
  xbin = hhpt.FindBin(datapt[ii-1][0])
  hhpt.SetBinContent(xbin, datapt[ii-1][1])

outFile.Write()
outFile.Close()
