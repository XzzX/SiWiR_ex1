import numpy as np
import matplotlib.pyplot as pl

dim1BLAS, dim2BLAS, dim3BLAS, tBLAS, L2BandwidthBLAS, L2RequestRateBLAS, L2MissRateBLAS, L2MissRatioBLAS, MFlopssBLAS, AVXMFlopssBLAS, PackedMUOPSsBLAS, ScalarMUOPSsBLAS  = np.loadtxt("data/perfBLAS.txt", skiprows=1).transpose()
dim1Naive, dim2Naive, dim3Naive, tNaive, L2BandwidthNaive, L2RequestRateNaive, L2MissRateNaive, L2MissRatioNaive, MFlopssNaive, AVXMFlopssNaive, PackedMUOPSsNaive, ScalarMUOPSsNaive  = np.loadtxt("data/perfNaive.txt", skiprows=1).transpose()
dim1Opt1, dim2Opt1, dim3Opt1, tOpt1, L2BandwidthOpt1, L2RequestRateOpt1, L2MissRateOpt1, L2MissRatioOpt1, MFlopssOpt1, AVXMFlopssOpt1, PackedMUOPSsOpt1, ScalarMUOPSsOpt1  = np.loadtxt("data/perfOpt1.txt", skiprows=1).transpose()
dim1Opt2, dim2Opt2, dim3Opt2, tOpt2, L2BandwidthOpt2, L2RequestRateOpt2, L2MissRateOpt2, L2MissRatioOpt2, MFlopssOpt2, AVXMFlopssOpt2, PackedMUOPSsOpt2, ScalarMUOPSsOpt2  = np.loadtxt("data/perfOpt2.txt", skiprows=1).transpose()
dim1Opt3, dim2Opt3, dim3Opt3, tOpt3, L2BandwidthOpt3, L2RequestRateOpt3, L2MissRateOpt3, L2MissRatioOpt3, MFlopssOpt3, AVXMFlopssOpt3, PackedMUOPSsOpt3, ScalarMUOPSsOpt3  = np.loadtxt("data/perfOpt3.txt", skiprows=1).transpose()
dim1Opt4, dim2Opt4, dim3Opt4, tOpt4, L2BandwidthOpt4, L2RequestRateOpt4, L2MissRateOpt4, L2MissRatioOpt4, MFlopssOpt4, AVXMFlopssOpt4, PackedMUOPSsOpt4, ScalarMUOPSsOpt4  = np.loadtxt("data/perfOpt4.txt", skiprows=1).transpose()
dim1Opt5, dim2Opt5, dim3Opt5, tOpt5, L2BandwidthOpt5, L2RequestRateOpt5, L2MissRateOpt5, L2MissRatioOpt5, MFlopssOpt5, AVXMFlopssOpt5, PackedMUOPSsOpt5, ScalarMUOPSsOpt5  = np.loadtxt("data/perfOpt5.txt", skiprows=1).transpose()

ind = np.arange(5, 12, 1)
width = 0.1

fig, ax = pl.subplots()
rectsBlasTime = ax.bar(ind-3*width, tBLAS, width, color='r')
rectsNaiveTime = ax.bar(ind-2*width, tNaive, width, color='y')
rectsOpt1Time = ax.bar(ind-width, tOpt1, width, color='b')
rectsOpt2Time = ax.bar(ind, tOpt2, width, color='g')
rectsOpt3Time = ax.bar(ind+width, tOpt3, width, color='k')
rectsOpt4Time = ax.bar(ind+2*width, tOpt4, width, color='w')
rectsOpt5Time = ax.bar(ind+3*width, tOpt5, width, color='c')

ax.set_title('Time Performance')
ax.set_ylabel('log(time [s])')
pl.yscale('log', nonposy='clip')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasTime[0], rectsNaiveTime[0], rectsOpt1Time[0], rectsOpt2Time[0], rectsOpt3Time[0], rectsOpt4Time[0], rectsOpt5Time[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/TimePerformance.pdf")

fig, ax = pl.subplots()
rectsBlasBandwidth = ax.bar(ind, L2BandwidthBLAS, width, color='r')
rectsNaiveBandwidth = ax.bar(ind+width, L2BandwidthNaive, width, color='y')
rectsOpt1Bandwidth = ax.bar(ind+2*width, L2BandwidthOpt1, width, color='b')
rectsOpt2Bandwidth = ax.bar(ind+3*width, L2BandwidthOpt2, width, color='g')
rectsOpt3Bandwidth = ax.bar(ind+4*width, L2BandwidthOpt3, width, color='k')
rectsOpt4Bandwidth = ax.bar(ind+5*width, L2BandwidthOpt4, width, color='w')
rectsOpt5Bandwidth = ax.bar(ind+6*width, L2BandwidthOpt5, width, color='c')
ax.set_title('L2 Bandwidth')
ax.set_ylabel('L2 bandwidth [MBytes/s]')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasBandwidth[0], rectsNaiveBandwidth[0], rectsOpt1Bandwidth[0], rectsOpt2Bandwidth[0], rectsOpt3Bandwidth[0], rectsOpt4Bandwidth[0], rectsOpt5Bandwidth[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/L2Bandwidth.pdf")

fig, ax = pl.subplots()
rectsBlasRequestRate = ax.bar(ind, L2RequestRateBLAS, width, color='r')
rectsNaiveRequestRate = ax.bar(ind+width, L2RequestRateNaive, width, color='y')
rectsOpt1RequestRate = ax.bar(ind+2*width, L2RequestRateOpt1, width, color='b')
rectsOpt2RequestRate = ax.bar(ind+3*width, L2RequestRateOpt2, width, color='g')
rectsOpt3RequestRate = ax.bar(ind+4*width, L2RequestRateOpt3, width, color='k')
rectsOpt4RequestRate = ax.bar(ind+5*width, L2RequestRateOpt4, width, color='w')
rectsOpt5RequestRate = ax.bar(ind+6*width, L2RequestRateOpt5, width, color='c')
ax.set_title('L2 Request Rate')
ax.set_ylabel('L2 request rate')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasRequestRate[0], rectsNaiveRequestRate[0], rectsOpt1RequestRate[0], rectsOpt2RequestRate[0], rectsOpt3RequestRate[0], rectsOpt4RequestRate[0], rectsOpt5RequestRate[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/L2RequestRate.pdf")

fig, ax = pl.subplots()
rectsBlasMissRate = ax.bar(ind, L2MissRateBLAS, width, color='r')
rectsNaiveMissRate = ax.bar(ind+width, L2MissRateNaive, width, color='y')
rectsOpt1MissRate = ax.bar(ind+2*width, L2MissRateOpt1, width, color='b')
rectsOpt2MissRate = ax.bar(ind+3*width, L2MissRateOpt2, width, color='g')
rectsOpt3MissRate = ax.bar(ind+4*width, L2MissRateOpt3, width, color='k')
rectsOpt4MissRate = ax.bar(ind+5*width, L2MissRateOpt4, width, color='w')
rectsOpt5MissRate = ax.bar(ind+6*width, L2MissRateOpt5, width, color='c')
ax.set_title('L2 Miss Rate')
ax.set_ylabel('L2 Miss rate')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasMissRate[0], rectsNaiveMissRate[0], rectsOpt1MissRate[0], rectsOpt2MissRate[0], rectsOpt3MissRate[0], rectsOpt4MissRate[0], rectsOpt5MissRate[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/L2MissRate.pdf")

fig, ax = pl.subplots()
rectsBlasMissRatio = ax.bar(ind, L2MissRatioBLAS, width, color='r')
rectsNaiveMissRatio = ax.bar(ind+width, L2MissRatioNaive, width, color='y')
rectsOpt1MissRatio = ax.bar(ind+2*width, L2MissRatioOpt1, width, color='b')
rectsOpt2MissRatio = ax.bar(ind+3*width, L2MissRatioOpt2, width, color='g')
rectsOpt3MissRatio = ax.bar(ind+4*width, L2MissRatioOpt3, width, color='k')
rectsOpt4MissRatio = ax.bar(ind+5*width, L2MissRatioOpt4, width, color='w')
rectsOpt5MissRatio = ax.bar(ind+6*width, L2MissRatioOpt5, width, color='c')
ax.set_title('L2 Miss Ratio')
ax.set_ylabel('L2 Miss Ratio')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasMissRatio[0], rectsNaiveMissRatio[0], rectsOpt1MissRatio[0], rectsOpt2MissRatio[0], rectsOpt3MissRatio[0], rectsOpt4MissRatio[0], rectsOpt5MissRatio[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/L2MissRatio.pdf")

fig, ax = pl.subplots()
rectsBlasMFlopss = ax.bar(ind, MFlopssBLAS, width, color='r')
rectsNaiveMFlopss = ax.bar(ind+width, MFlopssNaive, width, color='y')
rectsOpt1MFlopss = ax.bar(ind+2*width, MFlopssOpt1, width, color='b')
rectsOpt2MFlopss = ax.bar(ind+3*width, MFlopssOpt2, width, color='g')
rectsOpt3MFlopss = ax.bar(ind+4*width, MFlopssOpt3, width, color='k')
rectsOpt4MFlopss = ax.bar(ind+5*width, MFlopssOpt4, width, color='w')
rectsOpt5MFlopss = ax.bar(ind+6*width, MFlopssOpt5, width, color='c')
ax.set_title('MFlops/s')
ax.set_ylabel('MFlops/s')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasMFlopss[0], rectsNaiveMFlopss[0], rectsOpt1MFlopss[0], rectsOpt2MFlopss[0], rectsOpt3MFlopss[0], rectsOpt4MFlopss[0], rectsOpt5MFlopss[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/MFlopss.pdf")

fig, ax = pl.subplots()
rectsBlasAVXMFlopss = ax.bar(ind, AVXMFlopssBLAS, width, color='r')
rectsNaiveAVXMFlopss = ax.bar(ind+width, AVXMFlopssNaive, width, color='y')
rectsOpt1AVXMFlopss = ax.bar(ind+2*width, AVXMFlopssOpt1, width, color='b')
rectsOpt2AVXMFlopss = ax.bar(ind+3*width, AVXMFlopssOpt2, width, color='g')
rectsOpt3AVXMFlopss = ax.bar(ind+4*width, AVXMFlopssOpt3, width, color='k')
rectsOpt4AVXMFlopss = ax.bar(ind+5*width, AVXMFlopssOpt4, width, color='w')
rectsOpt5AVXMFlopss = ax.bar(ind+6*width, AVXMFlopssOpt5, width, color='c')
ax.set_title('AVXMFlops/s')
ax.set_ylabel('AVXMFlops/s')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasAVXMFlopss[0], rectsNaiveAVXMFlopss[0], rectsOpt1AVXMFlopss[0], rectsOpt2AVXMFlopss[0], rectsOpt3AVXMFlopss[0], rectsOpt4AVXMFlopss[0], rectsOpt5AVXMFlopss[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/AVXMFlopss.pdf")

fig, ax = pl.subplots()
rectsBlasPackedMUOPSs = ax.bar(ind, PackedMUOPSsBLAS, width, color='r')
rectsNaivePackedMUOPSs = ax.bar(ind+width, PackedMUOPSsNaive, width, color='y')
rectsOpt1PackedMUOPSs = ax.bar(ind+2*width, PackedMUOPSsOpt1, width, color='b')
rectsOpt2PackedMUOPSs = ax.bar(ind+3*width, PackedMUOPSsOpt2, width, color='g')
rectsOpt3PackedMUOPSs = ax.bar(ind+4*width, PackedMUOPSsOpt3, width, color='k')
rectsOpt4PackedMUOPSs = ax.bar(ind+5*width, PackedMUOPSsOpt4, width, color='w')
rectsOpt5PackedMUOPSs = ax.bar(ind+6*width, PackedMUOPSsOpt5, width, color='c')
ax.set_title('PackedMUOPS/s')
ax.set_ylabel('PackedMUOPS/s')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasPackedMUOPSs[0], rectsNaivePackedMUOPSs[0], rectsOpt1PackedMUOPSs[0], rectsOpt2PackedMUOPSs[0], rectsOpt3PackedMUOPSs[0], rectsOpt4PackedMUOPSs[0], rectsOpt5PackedMUOPSs[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/PackedMUOPSs.pdf")

fig, ax = pl.subplots()
rectsBlasScalarMUOPSs = ax.bar(ind, ScalarMUOPSsBLAS, width, color='r')
rectsNaiveScalarMUOPSs = ax.bar(ind+width, ScalarMUOPSsNaive, width, color='y')
rectsOpt1ScalarMUOPSs = ax.bar(ind+2*width, ScalarMUOPSsOpt1, width, color='b')
rectsOpt2ScalarMUOPSs = ax.bar(ind+3*width, ScalarMUOPSsOpt2, width, color='g')
rectsOpt3ScalarMUOPSs = ax.bar(ind+4*width, ScalarMUOPSsOpt3, width, color='k')
rectsOpt4ScalarMUOPSs = ax.bar(ind+5*width, ScalarMUOPSsOpt4, width, color='w')
rectsOpt5ScalarMUOPSs = ax.bar(ind+6*width, ScalarMUOPSsOpt5, width, color='c')
ax.set_title('ScalarMUOPS/s')
ax.set_ylabel('ScalarMUOPS/s')
ax.set_xticks(ind+width)
ax.set_xticklabels(('5', '6', '7', '8', '9', '10', '11'))
ax.set_xlabel('log2(matrix size)')
ax.legend((rectsBlasScalarMUOPSs[0], rectsNaiveScalarMUOPSs[0], rectsOpt1ScalarMUOPSs[0], rectsOpt2ScalarMUOPSs[0], rectsOpt3ScalarMUOPSs[0], rectsOpt4ScalarMUOPSs[0], rectsOpt5ScalarMUOPSs[0]), ('BLAS', 'Naive', 'Opt1', 'Opt2', 'Opt3', 'Opt4', 'Opt5'), loc=2)
pl.savefig("data/PackedMUOPSs.pdf")