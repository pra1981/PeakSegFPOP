four.out: PeakSegFPOP four.bedGraph
	./PeakSegFPOP four.bedGraph 0.1 | tee four.out
test: test.R
	R --no-save < $<
three.out: PeakSegFPOP three.bedGraph
	./PeakSegFPOP three.bedGraph 0.1 | tee three.out
PeakSegFPOP: *.cpp *.h 
	g++ -std=c++0x -o PeakSegFPOP PeakSegFPOPLog.cpp funPieceListLog.cpp -ldb_stl
hg19_problems.bed: gap2problems.R hg19_gap.bed hg19_chromInfo.txt 
	Rscript gap2problems.R hg19_gap.bed hg19_chromInfo.txt hg19_problems.bed
