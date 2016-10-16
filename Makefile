four.out: PeakSegFPOP four.bedGraph
	./PeakSegFPOP four.bedGraph 0.1 | tee four.out
test_run.out: test_run.R test_cases.R PeakSegFPOP
	cp PeakSegFPOP ~/bin
	R --vanilla < test_run.R
	touch test_run.out
three.out: PeakSegFPOP three.bedGraph
	./PeakSegFPOP three.bedGraph 0.1 | tee three.out
PeakSegFPOP: *.cpp *.h 
	g++ -std=c++0x -o PeakSegFPOP PeakSegFPOPLog.cpp funPieceListLog.cpp -ldb_stl
hg19_problems.bed: gap2problems.R hg19_gap.bed hg19_chromInfo.txt 
	Rscript gap2problems.R hg19_gap.bed hg19_chromInfo.txt hg19_problems.bed
hg19_problems_some.bed: hg19_problems.bed
	intersectBed -sorted -wa -u -a hg19_problems.bed -b labels/H3K36me3_AM_immune/McGill0322/coverage.bedGraph| tee hg19_problems_some.bed
